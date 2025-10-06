/*---------------------------------------------------------------------------*\
                        _   _ ____________ ___________
                       | | | ||  ___|  _  \_   _| ___ \     H ybrid
  ___  _ __   ___ _ __ | |_| || |_  | | | | | | | |_/ /     F ictitious
 / _ \| '_ \ / _ \ '_ \|  _  ||  _| | | | | | | | ___ \     D omain
| (_) | |_) |  __/ | | | | | || |   | |/ / _| |_| |_/ /     I mmersed
 \___/| .__/ \___|_| |_\_| |_/\_|   |___/  \___/\____/      B oundary
      | |
      |_|
-------------------------------------------------------------------------------
License

    openHFDIB is licensed under the GNU LESSER GENERAL PUBLIC LICENSE (LGPL).

    Everyone is permitted to copy and distribute verbatim copies of this license
    document, but changing it is not allowed.

    This version of the GNU Lesser General Public License incorporates the terms
    and conditions of version 3 of the GNU General Public License, supplemented
    by the additional permissions listed below.

    You should have received a copy of the GNU Lesser General Public License
    along with openHFDIB. If not, see <http://www.gnu.org/licenses/lgpl.html>.

InNamspace
    Foam

Contributors
    Federico Municchi (2018)
\*---------------------------------------------------------------------------*/
#include "immersedBody.H"
#include "fvMesh.H"
#include "polyMesh.H"
#include "fvCFD.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include <cmath>
#include <algorithm>

#include "interpolationCellPoint.H"
#include "interpolationCell.H"
#include "meshSearch.H"
#include "List.H"
#include "UPstream.H"
#include "transformGeometricField.H"
#include "triSurfaceSearch.H"
#include "triSurface.H"
#include "triSurfaceMesh.H"
#include "Tuple2.H"
#include "mathematicalConstants.H"
#include "pointIndexHit.H"



#define ORDER 2

using namespace Foam;

//---------------------------------------------------------------------------//
immersedBody::immersedBody
(
    word bodyName,
    const Foam::fvMesh& mesh,
    dictionary& HFDIBDict,
    dictionary& transportProperties
)
:
bodyName_(bodyName),
isFirstUpdate_(true),
immersedDict_(HFDIBDict.subDict(bodyName)),
mesh_(mesh),
transportProperties_(transportProperties),
M_(0.),
CoM_(vector::zero),
Axis_(vector::zero),
omega_(0.),
Vel_(vector::zero),
scaleEnabled_(false),   // <-- NEW
t0_(0),                 // <-- NEW
t1_(0),                 // <-- NEW
s0_(1),                 // <-- NEW
s1_(1),                 // <-- NEW
sPrev_(1),            // <— NEW
I_(symmTensor::zero),
bodySurfMesh_
(
    new triSurfaceMesh
    (
        IOobject
        (
	 fileName(immersedDict_.lookup("fileName")), //changed
            mesh_.time().constant(),
            "triSurface",
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    )
)
{

    //Process stl file
    Info<< "Read Immersed Boundary triSurface "
        << immersedDict_.lookup("fileName")
        << " for body " << bodyName << endl;

    bodySurfMesh_->writeStats(Info);
    Info << endl;

    if(immersedDict_.found("transform"))
    {
        dictionary transformDict = immersedDict_.subDict("transform");

        transformBody(transformDict);
    }

    //Return if declared static
    if(immersedDict_.found("staticBody") )
    {
        bodyOperation_=STATICBODY;
        Info << bodyName << " is static body." << endl;
    }
    else if(immersedDict_.found("transRotatingBody"))
    {
        bodyOperation_=TRANSROTATINGBODY;

        //Get basic quantities from dict
        const dictionary& tr = immersedDict_.subDict("transRotatingBody");

        // Axis_ = tr.lookup<vector>("axis");
        // CoM_  = tr.lookup<vector>("center");
        // omega_ = readScalar(tr.lookup("omega"));
        // Vel_  = tr.lookup<vector>("velocity");

	 // AFTER (works in v2412)
	 tr.lookup("axis")     >> Axis_;
	 tr.lookup("center")   >> CoM_;
	 omega_ = readScalar(tr.lookup("omega"));
	 tr.lookup("velocity") >> Vel_;

	// --- NEW: optional scaling over time ---

	// Optional: normalize axis to be safe
       Axis_ /= Foam::max(VSMALL, mag(Axis_));

       //	if (immersedDict_.subDict("transRotatingBody").found("scale"))
       //  {
       //   const dictionary sD =
       //     immersedDict_.subDict("transRotatingBody").subDict("scale");
       //
       //   t0_ = sD.lookupOrDefault<scalar>("t0", scalar(0));
       //   t1_ = sD.lookupOrDefault<scalar>("t1", scalar(0));
       //   s0_ = sD.lookupOrDefault<scalar>("s0", scalar(1));
       //  s1_ = sD.lookupOrDefault<scalar>("s1", scalar(1));
       //
	    //  scaleEnabled_ = true;
	    //	      }
       if (tr.found("scale"))
	 {
	   const dictionary sD = tr.subDict("scale");
	   t0_ = sD.lookupOrDefault<scalar>("t0", scalar(0));
	   t1_ = sD.lookupOrDefault<scalar>("t1", scalar(0));
	   s0_ = sD.lookupOrDefault<scalar>("s0", scalar(1));
	   s1_ = sD.lookupOrDefault<scalar>("s1", scalar(1));
	   scaleEnabled_ = true;
	   sPrev_ = 1.0;
	 }

        Info << bodyName << " has scripted trans rotational motion." << endl;
    }
    else if(immersedDict_.found("fluidCoupling"))
    {
        bodyOperation_=FLUIDCOUPLING;
        Info << bodyName << " is coupled with fluid phase." << endl;

    }
    else
    {
        Info << "No body operation was found for " << bodyName << endl
             << "Assuming static body.";
        bodyOperation_=STATICBODY;
    }

    // --- leafletRT: optional prescribed motion block -----------------
    if (immersedDict_.found("deformingBody"))
      {
	const dictionary& d = immersedDict_.subDict("deformingBody");
	const word t = d.lookupOrDefault<word>("type", "leafletRT");

	if (t == "leafletRT")
	  {
	    readLeafletRT(d);
	    if (autoR0_) computeAutoR0_();
	    Info<< bodyName << " has leafletRT prescribed motion." << nl;
	  }
	else if (t == "leafletPlane")
	  {
	    readLeafletPlane(d);
	    if (autoH0_) computeAutoH0_();
	    Info<< bodyName << " has leafletPlane prescribed motion." << nl;
	  }
	else if (t == "annulusRZ")
	  {
	    readAnnulusRZ(d);
	    annulusEnabled_ = true;
	    Info<< bodyName << " has annulusRZ prescribed motion." << nl;
	  }
	else
	  {
	    FatalErrorInFunction
	      << "Type = " << t << " is unknown. Available options are "
	      << "'leafletRT', 'leafletPlane', and 'annulusRZ'" << exit(FatalError);
	  }
      }

    if (immersedDict_.found("deformingBody_mitral"))
    {
        const dictionary& d = immersedDict_.subDict("deformingBody_mitral");
        readMitralSliceAxis(d);
        mitralEnabled_ = true;
        Info<< bodyName_ << " : deformingBody_mitral (slice-axis) enabled" << nl;
    }
}
//---------------------------------------------------------------------------//
immersedBody::~immersedBody()
{
    bodySurfMesh_.clear();
}
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
// leafletRT helpers
// piecewise-linear interpolation
scalar immersedBody::interp(const List<scalar>& x, const List<scalar>& y, scalar xv) const
{
    if (x.empty()) return 0.0;
    if (xv <= x.first()) return y.first();
    if (xv >= x.last())  return y.last();
    label i = 1;
    for (; i < x.size(); ++i) if (xv < x[i]) break;
    const scalar x0 = x[i-1], x1 = x[i];
    const scalar y0 = y[i-1], y1 = y[i];
    const scalar a  = (xv - x0)/(x1 - x0 + VSMALL);
    return (1.0 - a)*y0 + a*y1;
}

// compute R0_ from STL points as max radius in plane ⟂ leafletAxis_ about leafletCenter_
void immersedBody::computeAutoR0_()
{
    if (!autoR0_) return;

    // const pointField pts(bodySurfMesh_->points());

    tmp<pointField> tPts = bodySurfMesh_->points();
    const pointField& pts = tPts();

    vector n = leafletAxis_;
    n /= (mag(n) + VSMALL);

    scalar RmaxFound = 0.0;
    forAll(pts, i)
    {
        vector rp = vector(pts[i] - leafletCenter_);
        scalar zc = rp & n;
        vector inPlane = rp - zc*n;
        RmaxFound = Foam::max(RmaxFound, mag(inPlane));
    }

    R0_  = RmaxFound;
    if (Rmax_ <= VSMALL) Rmax_ = R0_;  // default: clamp at rim
    Info<< "leafletRT: autoR0 computed from STL for body " << bodyName_
        << " => R0 = " << R0_ << " m" << nl;
}

// read deformingBody{ type leafletRT; ... } block
void immersedBody::readLeafletRT(const dictionary& d)
{
    leafletEnabled_ = true;

    leafletAxis_   = vector(d.lookup("axis"));
    leafletAxis_  /= (mag(leafletAxis_) + VSMALL);

    leafletCenter_ = point(d.lookup("center"));

    autoR0_ = d.lookupOrDefault<bool>("autoR0", false);

    if (!autoR0_)
    {
        if (!d.found("R0"))
        {
            FatalErrorInFunction
                << "For leafletRT: either set 'autoR0 true;' or supply 'R0 <value>;'"
                << exit(FatalError);
        }
        R0_ = readScalar(d.lookup("R0"));
    }

    Rmin_ = d.lookupOrDefault<scalar>("Rmin", 0.0);
    Rmax_ = d.lookupOrDefault<scalar>("Rmax", 0.0);   // 0 => use R0 at runtime
    clampOutside_ = d.lookupOrDefault<bool>("clampOutside", true);

    // radialLaw { table ( (rhat f) ... ) }
    const dictionary& rl = d.subDict("radialLaw");
    {
        const List<List<scalar> > tab(rl.lookup("table"));
        rHat_.setSize(tab.size());
        fR_.setSize(tab.size());
        forAll(tab, i){ rHat_[i] = tab[i][0]; fR_[i] = tab[i][1]; }
    }

    // amplitude & offset tables
    auto readPW = [&](const word& name, List<scalar>& tx, List<scalar>& ty)
    {
        const dictionary& s = d.subDict(name);
        const List<List<scalar> > tab(s.lookup("table"));
        tx.setSize(tab.size()); ty.setSize(tab.size());
        forAll(tab, i){ tx[i]=tab[i][0]; ty[i]=tab[i][1]; }
    };
    readPW("amplitude", tA_, A_);
    readPW("offset",    tB_, B_);

    Info<< "leafletRT: prescribed motion enabled for body " << bodyName_ << nl;
}

// apply Δz(r,t) = A(t)*f(R/R0) + B(t) along leafletAxis_
void immersedBody::applyLeafletRT()
{
    if (!leafletEnabled_) return;

    // ensure radius is known
    if (autoR0_ && (R0_ <= SMALL)) computeAutoR0_();
    if (R0_ <= SMALL)
    {
        FatalErrorInFunction
            << "leafletRT: R0_ is not set. Did you forget 'R0' or 'autoR0 true;'?"
            << exit(FatalError);
    }

    const scalar tNow = mesh_.time().value();
    const scalar A    = AofT(tNow);
    const scalar B    = BofT(tNow);

    // effective outer clamp
    const scalar Rcap = (Rmax_ > VSMALL ? Rmax_ : R0_);

    // pointField bodyPoints (bodySurfMesh_->points());

    tmp<pointField> tBody = bodySurfMesh_->points();
    pointField bodyPoints(tBody());

    vector n = leafletAxis_;
    n /= (mag(n) + VSMALL);

    forAll(bodyPoints, p)
    {
        vector rp  = vector(bodyPoints[p] - leafletCenter_);
        scalar zc  = rp & n;
        vector inP = rp - zc*n;
        const scalar R = mag(inP);

        if (clampOutside_)
        {
            if (R > Rcap + SMALL) continue;                 // freeze beyond rim
            if (R < Foam::max(0.0, Rmin_) - SMALL) continue;      // optional inner clamp
        }

        const scalar dz = A * fOfR(R) + B;
        bodyPoints[p] += dz * n;
    }

    bodySurfMesh_->movePoints(bodyPoints);
}
//------------------------------------------------------------------//

void immersedBody::readLeafletPlane(const dictionary& d)
{
    planeEnabled_ = true;

    // Motion direction (reuse leafletAxis_ key name "axis")
    leafletAxis_ = vector(d.lookup("axis"));
    leafletAxis_ /= (mag(leafletAxis_) + VSMALL);

    // Hinge plane: normal + a point on the plane
    planeNormal_ = vector(d.lookup("planeNormal"));
    planeNormal_ /= (mag(planeNormal_) + VSMALL);
    planePoint_  = point(d.lookup("planePoint"));

    // Auto span or fixed H0
    autoH0_ = d.lookupOrDefault<bool>("autoH0", false);
    if (!autoH0_)
    {
        if (!d.found("H0"))
        {
            FatalErrorInFunction
                << "leafletPlane: need 'H0 <value>;' or 'autoH0 true;'"
                << exit(FatalError);
        }
        H0_ = readScalar(d.lookup("H0"));
    }

    // Clamp handling (shared with leafletRT)
    clampOutside_ = d.lookupOrDefault<bool>("clampOutside", true);

    // Span law g(d̂) in [0..1] — reuse rHat_/fR_ storage
    const dictionary& sl = d.subDict("spanLaw");
    {
        const List<List<scalar> > tab(sl.lookup("table"));
        rHat_.setSize(tab.size());
        fR_.setSize(tab.size());
        forAll(tab, i){ rHat_[i] = tab[i][0]; fR_[i] = tab[i][1]; }
    }

    // Time tables — reuse tA_/A_, tB_/B_
    auto readPW = [&](const word& name, List<scalar>& tx, List<scalar>& ty)
    {
        const dictionary& s = d.subDict(name);
        const List<List<scalar> > tab(s.lookup("table"));
        tx.setSize(tab.size()); ty.setSize(tab.size());
        forAll(tab, i){ tx[i]=tab[i][0]; ty[i]=tab[i][1]; }
    };
    readPW("amplitude", tA_, A_);
    readPW("offset",    tB_, B_);
}

//-----------------------------------------------------------------------//

void immersedBody::computeAutoH0_()
{
    if (!planeEnabled_) return;

    tmp<pointField> tPts = bodySurfMesh_->points();
    const pointField& pts = tPts();

    scalar maxAbs = 0.0;
    forAll(pts, i)
    {
        const scalar d = (vector(pts[i] - planePoint_) & planeNormal_);
        maxAbs = max(maxAbs, mag(d));
    }
    H0_ = max(SMALL, maxAbs);
}
//---------------------------------------------------------------------//
void immersedBody::applyLeafletPlane()
{
  // Info<< __FILE__ << " " << __LINE__ << endl;
    if (!planeEnabled_) return;
  // Info<< __FILE__ << " " << __LINE__ << endl;

    if (autoH0_ && (H0_ <= SMALL)) computeAutoH0_();
    if (H0_ <= SMALL)
    {
        FatalErrorInFunction
            << "leafletPlane: H0_ not set. Need 'H0' or 'autoH0 true;'."
            << exit(FatalError);
    }

    const scalar tNow = mesh_.time().value();
    const scalar A    = AofT(tNow);
    const scalar B    = BofT(tNow);

    tmp<pointField> tBody = bodySurfMesh_->points();
    pointField bodyPoints(tBody());

    // Normalised directions
    vector ez = leafletAxis_;  ez /= (mag(ez) + VSMALL);
    vector ny = planeNormal_;  ny /= (mag(ny) + VSMALL);

    // Make planeNormal sign irrelevant: detect which side the leaflet lies on
    scalar sumSigned = 0.0;
    forAll(bodyPoints, p)
    {
        sumSigned += (vector(bodyPoints[p] - planePoint_) & ny);
    }
    const scalar sgn = (sumSigned < 0.0 ? -1.0 : 1.0); // expect negative if leaflet is below plane for ny=(0 1 0)

    // --- DEBUG: deformation stats (begin) ---
    scalar minDisp =  GREAT;
    scalar maxDisp = -GREAT;
    label  movedCt = 0;
    // --- DEBUG: deformation stats (end) ---

    forAll(bodyPoints, p)
    {
        const scalar dSigned = (vector(bodyPoints[p] - planePoint_) & ny);
        const scalar dEff    = sgn * dSigned;   // >=0 on the leaflet side, 0 at the plane

        if (clampOutside_)
        {
            if (dEff < -SMALL)           continue;    // opposite side of plane
            if (dEff > H0_ + SMALL)      continue;    // beyond tip span
        }

        scalar dhat = dEff / (H0_ + VSMALL);          // 0 at hinge, 1 at tip
        dhat = max(0.0, min(1.0, dhat));

        const scalar g    = interp(rHat_, fR_, dhat); // span gain g(d̂)
        const scalar disp = A*g + B;                  // metres along ez

	// --- DEBUG: track stats ---
	minDisp = Foam::min(minDisp, disp);
	maxDisp = Foam::max(maxDisp, disp);
	++movedCt;
	// --- DEBUG end ---

        bodyPoints[p] += disp * ez;
    }

    if ( (mesh_.time().timeIndex() % 50) == 0 )
      {
	Info<< "leafletPlane: t=" << mesh_.time().value()
	    << " A(t)=" << A << " B(t)=" << B
	    << " movedPts=" << movedCt
	    << " disp[min,max]=[" << minDisp << ", " << maxDisp << "] m"
	    << " (H0=" << H0_ << ", axis=" << leafletAxis_ << ")"
	    << nl;
      }

    bodySurfMesh_->movePoints(bodyPoints);
  // Info<< __FILE__ << " " << __LINE__ << endl;
}
//---------------------------------------------------------------------//
// Robust TOP-annulus plane fit from the body surface.
// - Clusters faces by near-equal normals
// - Picks the flattest, large-area cluster (a cap)
// - If the opposite cap is present, chooses the TOP one deterministically
//   by signed distance along the normal from the global surface centroid.
// Outputs:
//   nHat     : unit normal of the TOP cap (annulus plane)
//   centroid : area-weighted centroid of the TOP cap
void immersedBody::autoFitPlaneFromSurface(vector& nHat, point& centroid) const
{
    // --- Input geometry ---
    const triSurface& ts = bodySurfMesh_->surface();
    tmp<pointField> tP = bodySurfMesh_->points();
    const pointField& P = tP();

    if (ts.size() == 0 || P.size() == 0)
    {
        Warning<< "mitral: empty triSurface or points; using defaults." << nl;
        nHat = vector(0,0,1);
        centroid = point::zero;
        return;
    }

    // --- Collect per-face data ---
    struct FaceRec { vector n; scalar A; point c; label f; };
    List<FaceRec> rec(ts.size());

    forAll(ts, fi)
    {
        const labelledTri& t = ts[fi];
        const point& pa = P[t[0]];
        const point& pb = P[t[1]];
        const point& pc = P[t[2]];

        const triPointRef tri(pa, pb, pc);
        const vector AN = tri.areaNormal();          // 2*area * unitNormal
        const scalar area = 0.5*mag(AN);
        vector n = vector::zero;
        if (area > SMALL) n = AN / (2.0*area);       // unit normal

        const point  c  = (pa + pb + pc)/3.0;

        rec[fi] = { n, area, c, fi };
    }

    // --- Cluster faces by near-equal normals ---
    const scalar cosTol = Foam::cos(10.0 * Foam::constant::mathematical::pi/180.0);
    // tighten to 6° if your STL is very clean
    List<vector>                 clusterN;     // running mean normal
    List<scalar>                 clusterA;     // total area
    List<DynamicList<label>>     clusterFaces; // member face indices

    forAll(rec, i)
    {
        if (rec[i].A <= SMALL || mag(rec[i].n) <= SMALL) continue;

        label hit = -1;
        forAll(clusterN, k)
        {
            if ((rec[i].n & clusterN[k]) > cosTol) { hit = k; break; }
        }

        if (hit < 0)
        {
            clusterN.append(rec[i].n);
            clusterA.append(rec[i].A);
            clusterFaces.append(DynamicList<label>(8));
            clusterFaces.last().append(i);
        }
        else
        {
            const scalar Aold = clusterA[hit];
            clusterA[hit] += rec[i].A;
            clusterN[hit]  = (Aold*clusterN[hit] + rec[i].A*rec[i].n) / clusterA[hit];
            clusterFaces[hit].append(i);
        }
    }

    if (clusterN.empty())
    {
        Warning<< "mitral: no valid face clusters; using defaults." << nl;
        nHat = vector(0,0,1);
        centroid = point::zero;
        return;
    }

    // --- Score clusters: prefer tight spread + large area (cap-like) ---
    label  best = -1;
    scalar bestScore = -GREAT;

    forAll(clusterN, k)
    {
        scalar cosMin = 1.0;
        for (label j : clusterFaces[k])
            cosMin = min(cosMin, clusterN[k] & rec[j].n);

        const scalar score = clusterA[k] * max(0.0, cosMin - cosTol);
        if (score > bestScore) { bestScore = score; best = k; }
    }

    if (best < 0)
    {
        Warning<< "mitral: unable to select a cap cluster; using defaults." << nl;
        nHat = vector(0,0,1);
        centroid = point::zero;
        return;
    }

    // --- Compute global area-weighted centroid (for top/bottom disambiguation) ---
    point Cglobal = point::zero;
    scalar Aglobal = 0.0;
    forAll(rec, i)
    {
        if (rec[i].A > SMALL)
        {
            Cglobal += rec[i].A * rec[i].c;
            Aglobal += rec[i].A;
        }
    }
    if (Aglobal > VSMALL) Cglobal /= Aglobal; else Cglobal = point::zero;

    // --- If an opposite cap exists, pick the TOP one deterministically ---
    // Opposite: normal ~ 180° from best
    const scalar cosOpp = Foam::cos(8.0  * Foam::constant::mathematical::pi/180.0);
    label bestOpp = -1;
    forAll(clusterN, k)
    {
        if (k == best) continue;
        if ((clusterN[k] & clusterN[best]) < -cosOpp) { bestOpp = k; break; }
    }

    // Build area-weighted centroids for the candidate(s)
    auto clusterCentroid = [&](label k)->point
    {
        point csum = point::zero; scalar At = 0.0;
        for (label j : clusterFaces[k]) { csum += rec[j].A * rec[j].c; At += rec[j].A; }
        return (At > VSMALL ? csum/At : rec[clusterFaces[k].first()].c);
    };

    label chosen = best;
    if (bestOpp != -1)
    {
        const point cBest    = clusterCentroid(best);
        const point cBestOpp = clusterCentroid(bestOpp);

        // Signed height along the best normal, relative to global centroid
        const scalar sBest    = vector(cBest    - Cglobal) & clusterN[best];
        const scalar sBestOpp = vector(cBestOpp - Cglobal) & clusterN[best];

        // We want the cap that sits "higher" along the (best) normal to be the TOP.
        // If opposite cap is higher along this normal, choose it as TOP.
        if (sBestOpp > sBest) chosen = bestOpp;
    }

    // --- Final outputs: area-weighted centroid & unit normal of the TOP cap ---
    centroid = clusterCentroid(chosen);
    const vector nChosen = clusterN[chosen];
    const scalar nmag = mag(nChosen);
    nHat = (nmag > VSMALL ? nChosen / nmag : vector(0,0,1));
}

// Use MitralValve.stl to decide the sign of the normal (axis direction)
void immersedBody::maybeFlipNormalUsingValve(vector& nHat, const point& centroid) const
{
    if (!m_flipNormalTowardValve_ || m_valveRefFile_.empty()) return;

    // Load reference STL just to get its centroid
    triSurface tsRef(m_valveRefFile_);
    if (!tsRef.size()) return;

    // Compute reference centroid
    point cRef = Foam::point::zero;
    scalar A   = 0.0;
    forAll(tsRef, fi)
    {
        const labelledTri& tri = tsRef[fi];
        const point& a = tsRef.points()[tri[0]];
        const point& b = tsRef.points()[tri[1]];
        const point& c = tsRef.points()[tri[2]];
        const vector n = Foam::triPointRef(a,b,c).areaNormal();
        const scalar dA = 0.5*mag(n);
        const point cTri( (a+b+c)/3.0 );
        cRef += dA * cTri;
        A    += dA;
    }
    if (A > VSMALL) cRef /= A; else return;

    // If valve centroid lies opposite to nHat, flip so that nHat points toward valve
    const vector v = vector(cRef - centroid);
    if ( (v & nHat) < 0 ) nHat = -nHat;
}


// Robustly find TOP (annulus) and BOTTOM (free-edge) caps on the INNER surface
// using the same normal-clustering approach as autoFitPlaneFromSurface.
// Returns Ct (top centroid), Cb (bottom centroid), and axisUnit = normalize(Cb - Ct).
void immersedBody::fitTopBottomCapsFromSurface(point& Ct, point& Cb, vector& axisUnit) const
{
    const triSurface& ts = bodySurfMesh_->surface();
    tmp<pointField> tP = bodySurfMesh_->points();
    const pointField& P = tP();

    // Collect face normals, areas, centroids (same as your existing method)
    struct FaceRec { vector n; scalar A; point c; label f; };
    List<FaceRec> rec(ts.size());

    forAll(ts, fi)
    {
        const labelledTri& t = ts[fi];
        const point& pa = P[t[0]];
        const point& pb = P[t[1]];
        const point& pc = P[t[2]];
        const triPointRef tri(pa, pb, pc);
        const vector AN = tri.areaNormal();
        const scalar area = 0.5*mag(AN);
        vector n = vector::zero;
        if (area > SMALL) n = AN / (2.0*area);
        const point  c  = (pa + pb + pc)/3.0;
        rec[fi] = { n, area, c, fi };
    }

    // Cluster by normals
    const scalar cosTol = Foam::cos(10.0 * Foam::constant::mathematical::pi/180.0);
    List<vector>                 clusterN;
    List<scalar>                 clusterA;
    List<DynamicList<label>>     clusterFaces;

    forAll(rec, i)
    {
        if (rec[i].A <= SMALL || mag(rec[i].n) <= SMALL) continue;
        label hit = -1;
        forAll(clusterN, k)
        {
            if ( (rec[i].n & clusterN[k]) > cosTol ) { hit = k; break; }
        }
        if (hit < 0)
        {
            clusterN.append(rec[i].n);
            clusterA.append(rec[i].A);
            clusterFaces.append(DynamicList<label>(8));
            clusterFaces.last().append(i);
        }
        else
        {
            const scalar Aold = clusterA[hit];
            clusterA[hit] += rec[i].A;
            clusterN[hit]  = (Aold*clusterN[hit] + rec[i].A*rec[i].n)/clusterA[hit];
            clusterFaces[hit].append(i);
        }
    }
    auto clusterCentroid = [&](label k)->point
    {
        point csum = point::zero; scalar At = 0.0;
        for (label j : clusterFaces[k]) { csum += rec[j].A * rec[j].c; At += rec[j].A; }
        return (At > VSMALL ? csum/At : rec[clusterFaces[k].first()].c);
    };

    // Choose the two most "cap-like" opposing clusters (largest area, near-opposite normals)
    label bestTop=-1, bestBot=-1;
    scalar bestAreaTop=-GREAT, bestAreaBot=-GREAT;

    forAll(clusterN, k)
    {
        // candidate normal
        const vector nk = (mag(clusterN[k])>VSMALL ? clusterN[k]/mag(clusterN[k]) : vector(0,0,1));
        // find the most opposite cluster
        label opp=-1; scalar cosMin=1.0;
        forAll(clusterN, j)
        {
            if (j==k) continue;
            const scalar c = nk & (mag(clusterN[j])>VSMALL ? clusterN[j]/mag(clusterN[j]) : vector(0,0,1));
            if (c < cosMin) { cosMin = c; opp = j; }
        }
        // keep the pair with the greatest combined area
        if (opp>=0 && cosMin < -0.95) // ~ >171°
        {
            const scalar Acomb = clusterA[k] + clusterA[opp];
            if (Acomb > bestAreaTop + bestAreaBot)
            {
                bestTop = k; bestBot = opp;
                bestAreaTop = clusterA[k]; bestAreaBot = clusterA[opp];
            }
        }
    }
    // Fallback: choose two largest-area clusters if no nearly-opposite pair found
    if (bestTop<0 || bestBot<0)
    {
        // sort by area
        List<label> idx(clusterA.size());
        forAll(idx,i) idx[i]=i;
        std::sort(idx.begin(), idx.end(), [&](label i, label j){return clusterA[i]>clusterA[j];});
        bestTop = (idx.size()>0 ? idx[0] : -1);
        bestBot = (idx.size()>1 ? idx[1] : -1);
    }

    if (bestTop<0 || bestBot<0)
    {
        WarningInFunction<< "fitTopBottomCapsFromSurface: could not find two cap clusters. Using default axis." << nl;
        Ct = point::zero; Cb = point(0,0,1);
        axisUnit = vector(0,0,1);
        return;
    }

    Ct = clusterCentroid(bestTop);
    Cb = clusterCentroid(bestBot);

    // Pick which is "top" by using your existing valve-orientation helper if available
    vector nTop = (mag(clusterN[bestTop])>VSMALL ? clusterN[bestTop]/mag(clusterN[bestTop]) : vector(0,0,1));
    maybeFlipNormalUsingValve(nTop, Ct);
    // Ensure axis points from TOP to BOTTOM consistently
    // After maybeFlip, if (Cb - Ct) · nTop < 0, swap (Ct,Cb)
    if ( (vector(Cb - Ct) & nTop) < 0 ) { point tmp=Ct; Ct=Cb; Cb=tmp; }

    axisUnit = vector(Cb - Ct);
    const scalar L = mag(axisUnit);
    axisUnit = (L>VSMALL ? axisUnit/L : vector(0,0,1));
}
//---------------------------------------------------------------------//
// ===== annulusRZ helpers & motion ==========================================
#include "mathematicalConstants.H"

// s(t) = sin^2(pi t / T)
scalar immersedBody::sOfT(const scalar t) const
{
    if (T_ <= VSMALL) return 0.0;
    const scalar a = Foam::sin(Foam::constant::mathematical::pi * t / T_);
    return a*a;
}

// u(z) = clamp((z - z0_)/zf, 0, 1)  -- z measured from base plane
scalar immersedBody::uOfZ(const scalar z) const
{
    if (zf_ <= VSMALL) return 1.0;
    const scalar u = (z - z0_) / zf_;
    return Foam::max(0.0, Foam::min(1.0, u));
}

// ri(z,t) using concave "ease-out" in z and sin^2 time easing:
// G(u) = 1 - (1 - u)^p,   u = clamp(z/zf, 0, 1)
scalar immersedBody::riOfZT(const scalar z, const scalar t) const
{
    const scalar s = sOfT(t);          // sin^2(pi t / T)
    const scalar u = uOfZ(z);          // 0..1
    const scalar G = 1.0 - Foam::pow(Foam::max(0.0, 1.0 - u), p_);
    return Foam::max(0.0, ri0_ * (1.0 - s*G));
}

scalar immersedBody::roOfZT(const scalar z, const scalar t) const
{
    return riOfZT(z,t) + thick_;
}
//----------------------------------------------------------------------//
// Parse deformingBody{ type annulusRZ; ... }
void immersedBody::readAnnulusRZ(const dictionary& d)
{
    annulusEnabled_ = true;

    // Axis & center
    annAxis_   = vector(d.lookup("axis"));
    annAxis_  /= (mag(annAxis_) + VSMALL);
    annCenter_ = point(d.lookup("center"));

    // Core parameters
    ri0_      = readScalar(d.lookup("ri0"));
    zf_       = readScalar(d.lookup("zf"));
    thick_    = d.lookupOrDefault<scalar>("thickness", 1.0);
    T_        = d.lookupOrDefault<scalar>("T", 1.0);
    tol_      = d.lookupOrDefault<scalar>("tol", 1e-6);
    wallTol_  = d.lookupOrDefault<scalar>("wallTol", 4e-4);
    p_        = d.lookupOrDefault<scalar>("p", 3.0);
    tolZ_     = d.lookupOrDefault<scalar>("tolZ", 1e-6);
    bandPolicy_ = d.lookupOrDefault<word>("bandPolicy", "baseFixed");

    const scalar ro0 = ri0_ + thick_;

    // --- REST geometry snapshot & classification ---
    tmp<pointField> tPts = bodySurfMesh_->points();
    X0_ = tPts(); // store rest points

    const label nPts = X0_.size();
    pointClass_.setSize(nPts, -1);
    alpha_.setSize(nPts, 0.0);
    isBase_.setSize(nPts, false);
    isTop_.setSize(nPts, false);

    scalar zMin = GREAT, zMax = -GREAT;

    forAll(X0_, i)
    {
        const vector rp = vector(X0_[i] - annCenter_);
        const scalar z  = rp & annAxis_;
        const vector inPlane = rp - z*annAxis_;
        const scalar R  = mag(inPlane);

        zMin = min(zMin, z);
        zMax = max(zMax, z);

        // Tag inner/outer/band by REST radius
        if (mag(R - ri0_) <= wallTol_)           { pointClass_[i] = 0; }
        else if (mag(R - ro0) <= wallTol_)       { pointClass_[i] = 2; }
        else if (R > ri0_ + SMALL && R < ro0 - SMALL)
        {
            pointClass_[i] = 1; // band
            alpha_[i] = clamp((R - ri0_)/(thick_ + VSMALL), 0.0, 1.0);
        }
        // else: untouched (-1)
    }

    z0_ = zMin;                  // base plane from REST STL
    L_  = zMax - zMin;

    // Tag base/top using tolZ_ around z0_ / z0_+L_
    forAll(X0_, i)
    {
        const scalar z = (vector(X0_[i] - annCenter_) & annAxis_);
        isBase_[i] = (mag(z - z0_)       <= tolZ_);
        isTop_[i]  = (mag(z - (z0_+L_))  <= tolZ_);
    }

    Info<< "annulusRZ: body " << bodyName_
        << " ri0=" << ri0_ << " thick=" << thick_
        << " zf=" << zf_ << " p=" << p_
        << " base z0=" << z0_ << " L=" << L_
        << " bandPolicy=" << bandPolicy_ << nl;
}

//---------------------------------------------------------------------------------//

// Apply ri(z,t)/ro(z,t) while keeping [ri0,ri0+thick] band fixed (if requested)
void immersedBody::applyAnnulusRZ()
{
    if (!annulusEnabled_) return;

    const scalar tNow = mesh_.time().value();

    // Start from REST points always
    pointField bodyPoints(X0_);

    label movedCt = 0;

    forAll(bodyPoints, i)
    {
        if (pointClass_[i] == -1) continue;     // untouched triangles

        // Decompose the REST point
        const vector rp0 = vector(X0_[i] - annCenter_);
        const scalar z0  = rp0 & annAxis_;                 // REST z
        const vector inPlane0 = rp0 - z0*annAxis_;
        const scalar R0  = mag(inPlane0);
        const vector rhat = (R0 > VSMALL ? inPlane0/R0 : vector::zero);

        // Base clamp: if on the base ring, do not move at all
        if (bandPolicy_ == "baseFixed" && isBase_[i]) { continue; }

        // Current target radii at this z0
        const scalar ri = riOfZT(z0, tNow);
        const scalar ro = ri + thick_;

        scalar Rnew = R0; // fallback

        if      (pointClass_[i] == 0) Rnew = ri;                // inner skin
        else if (pointClass_[i] == 2) Rnew = ro;                // outer skin
        else if (pointClass_[i] == 1) Rnew = ri + alpha_[i]*thick_; // band slides

        // Rebuild point: NO axial displacement
        bodyPoints[i] = annCenter_ + z0*annAxis_ + Rnew*rhat;
        ++movedCt;
    }

    if ( (mesh_.time().timeIndex() % 50) == 0 )
    {
        Info<< "annulusRZ(baseFixed): t=" << tNow << " movedPts=" << movedCt << nl;
    }

    bodySurfMesh_->movePoints(bodyPoints);
}

//----------------------------------------------------------------//
void immersedBody::readMitralSliceAxis(const dictionary& d)
{
    // Period & gains
    mPeriod_    = d.lookupOrDefault<scalar>("period",    1.0);
    mNearRamp_  = d.lookupOrDefault<scalar>("nearRamp",  0.004);
    mActiveLen_ = d.lookupOrDefault<scalar>("activeLen", 0.016);
    mFmax_      = d.lookupOrDefault<scalar>("Fmax",      0.995);
    mTimeLaw_   = d.lookupOrDefault<word>  ("timeLaw",   word("cos2"));
    tLeaf_      = d.lookupOrDefault<scalar>("leafThickness", 0.0); // 0 => single-skin
    mSpaceLaw_ = d.lookupOrDefault<word>("spaceLaw", word("smootherstep"));

    // Optional: help orient "top" cap using a reference STL
    m_flipNormalTowardValve_ = d.lookupOrDefault<Switch>("flipNormalTowardValve", false);
    m_valveRefFile_          = d.lookupOrDefault<fileName>("valveRefFile", fileName(""));

    // Optional hints
    const bool haveCt = d.found("topCenter");
    const bool haveCb = d.found("bottomCenter");
    if (haveCt) mCt_ = point(d.lookup("topCenter"));
    if (haveCb) mCb_ = point(d.lookup("bottomCenter"));

    // Snapshot REST points
    tmp<pointField> tPts = bodySurfMesh_->points();
    mX0_ = tPts();

    // If Ct/Cb not supplied, detect from inner surface (STL) automatically
    if (!haveCt || !haveCb)
    {
        fitTopBottomCapsFromSurface(mCt_, mCb_, mAxisUnit_);
    }
    else
    {
        mAxisUnit_ = vector(mCb_ - mCt_);
        const scalar L = mag(mAxisUnit_);
        if (L > VSMALL) mAxisUnit_ /= L; else mAxisUnit_ = vector(0,0,1);
    }

    mLen_ = mag(vector(mCb_ - mCt_));
    if (mLen_ <= VSMALL)
    {
        FatalErrorInFunction << "Mitral slice-axis: Ct and Cb coincide or invalid; cannot define axis." << exit(FatalError);
    }

    // Precompute per-vertex ξ and in-plane geometry relative to Ct+ξ â
    const label n = mX0_.size();
    mXi_.setSize(n, 0.0);
    mRperp0_.setSize(n, vector::zero);
    mR0_.setSize(n, 0.0);
    mAlpha_.setSize(n, 0.0);
    mNthru_.setSize(n, vector::zero); // left zero unless you add pairing

    forAll(mX0_, i)
    {
        const vector d0 = vector(mX0_[i] - mCt_);
        scalar xi = (d0 & mAxisUnit_);      // signed axial distance from Ct
        xi = clamp(xi, 0.0, mLen_);         // keep in [0,L]
        mXi_[i] = xi;

        // slice center C(ξ) and in-plane vector r_perp
        const vector Cxi = vector(mCt_) + xi*mAxisUnit_;
        vector r = vector(mX0_[i] - Cxi);
        r -= (r & mAxisUnit_) * mAxisUnit_; // project
        mRperp0_[i] = r;
        mR0_[i] = mag(r);
    }

    // If you can pair inner/outer vertices, fill mAlpha_ and mNthru_ here
    // For now, leave them default (single-surface); thickness enforcement is optional.
}

void immersedBody::applyMitralSliceAxis()
{
    if (!mitralEnabled_) return;

    const scalar tNow = mesh_.time().value();
    const scalar gainT = timeGain_(tNow); // 0..1, peak at T/2

    tmp<pointField> tPts = bodySurfMesh_->points();
    pointField pts = tPts(); // start from current (we'll rebuild from REST)

    // We REBUILD from REST every step to avoid drift:
    pts = mX0_;

    label movedCt = 0;
    scalar minF = GREAT, maxF = -GREAT;

    forAll(pts, i)
      {
	const scalar xi   = mXi_[i];                // 0..L
	const scalar gS   = spaceGain_(xi);         // spatial gate 0..1
	const scalar F    = clamp(mFmax_ * gS * gainT, 0.0, 0.999);  // local inward fraction

	const vector r0   = mRperp0_[i];            // in-plane vector from slice center at t=0
	const scalar R0   = mR0_[i];                // |r0|

	// per-slice center on axis at this xi (now we actually use it)
	const vector Cxi  = vector(mCt_) + xi*mAxisUnit_;

	vector Xnew = vector(mX0_[i]);              // start from REST by default
	if (R0 > VSMALL && F > SMALL)
	  {
	    // Express the point relative to Cxi: X = Cxi + r0 at t=0
	    // After inward contraction by factor F: r -> (1-F) * r0
	    Xnew = Cxi + (1.0 - F) * r0;
	    ++movedCt;
	    minF = Foam::min(minF, F);
	    maxF = Foam::max(maxF, F);
	  }

	pts[i] = point(Xnew);
      }

    if ( (mesh_.time().timeIndex() % 10) == 0 )
    {
        if (movedCt == 0) { minF = 0; maxF = 0; }
        Info<< "mitral(slice-axis): t=" << tNow
            << " timeGain=" << gainT
            << " movedPts=" << movedCt
            << " F[min,max]=[" << minF << "," << maxF << "]"
            << nl;
    }

    bodySurfMesh_->movePoints(pts);
}

//-----------------------------------------------------------------//
void immersedBody::transformBody(dictionary& transformDict)
{
    Info<< "Transforming immersed body " << bodyName_
        << " using dictionary" << endl;

    //pointField bodyPoints = bodySurfMesh_->points();

    tmp<pointField> tBody = bodySurfMesh_->points();
    pointField bodyPoints(tBody());

    vector transVec = transformDict.lookupOrDefault<vector>
    (
        "translate",
        vector::zero
    );

//- TODO: Implement!
//    vector rotVec = transformDict.lookupOrDefault<vector>
//    (
//        "rotate",
//        vector::zero
//    );
//
//    vector scaleVec = transformDict.lookupOrDefault<vector>
//    (
//        "scale",
//        vector::one
//    );

    forAll(bodyPoints,pointI)
    {
        bodyPoints[pointI] += transVec;
    }

     bodySurfMesh_->movePoints(bodyPoints);
}

//---------------------------------------------------------------------------//
//Update immersed body
void immersedBody::updateBodyField( volScalarField& body,
                                    volVectorField & f
                                 )
{
  // Info<< __FILE__ << " " << __LINE__ << endl;

    ensureStorageForCurrentMesh_();   // <— NEW


    if(isFirstUpdate_)
    {
  // Info<< __FILE__ << " " << __LINE__ << endl;
        createImmersedBody( body );
        isFirstUpdate_ = false;
    }
    else
    {
  // Info<< __FILE__ << " " << __LINE__ << endl;
        updateImmersedBody( body, f );
    }
  // Info<< __FILE__ << " " << __LINE__ << endl;

    body.correctBoundaryConditions();
}
//---------------------------------------------------------------------------//
//Create immersed body info
void immersedBody::createImmersedBody(volScalarField& body)
{
    // MUST be first
    ensureStorageForCurrentMesh_();

    triSurface        ibTemp(bodySurfMesh_());
    triSurfaceSearch  ibTriSurfSearch(ibTemp);

    const pointField& pp   = mesh_.points();
    const vector      span = 1.1*(mesh_.bounds().span() + vector::one*SMALL);

    intCells_.clear();
    surfCells_.clear();
    surfNorm_.clear();

    // --------- 1) Seed using vertex-inside + proximity fallback ----------
    forAll(mesh_.C(), cellI)
    {
        const labelList& vLbl = mesh_.cellPoints()[cellI];
        const pointField vPts(pp, vLbl);

        boolList vInside = ibTriSurfSearch.calcInside(vPts);

        // count "true" manually (no countIf in OF-2412)
        label nInside = 0;
        forAll(vInside, kk) if (vInside[kk]) ++nInside;

        // denom as plain label (avoid dimensioned<int>)
        const label denom = Foam::max(vPts.size(), label(1));
        body[cellI] = scalar(nInside) / scalar(denom);

        const bool anyInside = (nInside > 0);

        // proximity fallback when no vertex is inside (thin sheet case)
        bool nearSurf = false;
        if (!anyInside)
        {
            const point         c   = mesh_.C()[cellI];
            const pointIndexHit hit = ibTriSurfSearch.nearest(c, span);

            // cell length scale ~ cube-root of cell volume
            const scalar h = std::cbrt( Foam::max(mesh_.V()[cellI], VSMALL) );
            nearSurf = (hit.hit() && mag(hit.point() - c) <= 0.75*h);

            if (nearSurf)
            {
                body[cellI] = 0.5;            // neutral seed; MC will overwrite
            }
        }

        if      (body[cellI] >= 0.999)                        { intCells_.append(cellI); }
        else if (body[cellI] > VSMALL && body[cellI] < 0.999) { surfCells_.append(cellI); }
        else if (nearSurf)                                    { surfCells_.append(cellI); }
    }

    // --------- 2) Monte-Carlo refinement of surface cells ----------
    refineBody(body, ibTriSurfSearch, pp);    // uses refineMC_ (set in dict)

    // clamp numerically (keep as fractions; do NOT binarize)
    forAll(body, i) body[i] = Foam::max(0.0, Foam::min(1.0, body[i]));

    // --------- 3) Reclassify after refinement ----------
    intCells_.clear();
    const label initialCap = Foam::max(surfCells_.size(), label(128));
    DynamicList<label> newSurf(initialCap);

    forAll(body, cellI)
    {
        const scalar v = body[cellI];
        if      (v >= 0.999)                      { intCells_.append(cellI); }
        else if (v > VSMALL && v < 0.999)         { newSurf.append(cellI);   }
    }
    List<label> tmp = newSurf.shrink();
    surfCells_.transfer(tmp);

    // --------- 4) Downstream geometry for IB forcing ----------
    calculateInterpolationPoints(body, ibTriSurfSearch);
    calculateGeometricalProperties(body);

    Info<< "immersedBody: body min/max = "
        << gMin(body) << " / " << gMax(body) << nl;
}

//---------------------------------------------------------------------------//
//Update immersed body info
void immersedBody::updateImmersedBody
(
    volScalarField & body,
    volVectorField & f
)
{
  ensureStorageForCurrentMesh_();   // <— NEW


  //Check Operation to perform
  // If static AND no leaflet, we truly do nothing
  if (bodyOperation_==STATICBODY && !leafletEnabled_ && !planeEnabled_ && !annulusEnabled_ && !mitralEnabled_)
    {
      Info<< "leafletEnabled = " << leafletEnabled_ << endl;
      return;
    }

  // Info<< __FILE__ << " " << __LINE__ << endl;

  // Fluid-driven update (if requested)
  if (bodyOperation_==FLUIDCOUPLING)
    {
  // Info<< __FILE__ << " " << __LINE__ << endl;
      updateCoupling(body,f);
    }

  // Geometry-dependent fields must be recomputed after motion
  resetBody(body);

  // Apply motions:
  // 1) Prescribed leaflet law (if enabled)
  if (leafletEnabled_)
    {
  // Info<< __FILE__ << " " << __LINE__ << endl;
      applyLeafletRT();
    }

  // 1b) Plane-hinged leaflet (if enabled)
  // Info<< __FILE__ << " " << __LINE__ << endl;
  // Info<< "planeEnabled = " << planeEnabled_ << endl;
  if (planeEnabled_)
  {
  // Info<< __FILE__ << " " << __LINE__ << endl;
    applyLeafletPlane();
  }

   if (annulusEnabled_)
     {
     applyAnnulusRZ();
     }   // <<< add this line

   if (mitralEnabled_)
   {
       applyMitralSliceAxis();
   }

  // 2) Legacy trans-rotational motion (if requested)
  if (bodyOperation_==TRANSROTATINGBODY) { moveImmersedBody(); }

  //TODO: inefficient algorithm...
  // Info<< __FILE__ << " " << __LINE__ << endl;
  createImmersedBody(body);
  // Info<< __FILE__ << " " << __LINE__ << endl;
}

//---------------------------------------------------------------------------//
void immersedBody::updateCoupling
(
    volScalarField & body,
    volVectorField & f
)
{
     ensureStorageForCurrentMesh_();

    const uniformDimensionedVectorField g =
    mesh_.lookupObject<uniformDimensionedVectorField>("g");

    // const dimensionedScalar rhof(transportProperties_.lookup("rho"));
    const dimensionedScalar rhof("rho", transportProperties_);

    // dimensionedScalar rho_(immersedDict_.lookup("rho"));
    const dimensionedScalar rho_("rho", immersedDict_);

    vector F(vector::zero);
    vector T(vector::zero);

    //Calcualate viscous force and torque
    forAll(surfCells_,cell)
    {
        label cellI = surfCells_[cell];

        F +=  f[cellI] * mesh_.V()[cellI] * rhof.value();
        T +=  ( (mesh_.C()[cellI]-CoM_) ^ f[cellI] )
                 *mesh_.V()[cellI]* rhof.value();

    }

    reduce(F, sumOp<vector>());
    reduce(T, sumOp<vector>());

    //Update body linear velocity
    Vel_ += mesh_.time().deltaT().value()
            * (
                ((mag(M_) > VSMALL) ? (F/M_) : vector::zero)
                + (1.0-rhof.value()/rho_.value())*g.value()
              );

    //Update body angular velocity
    vector Omega_(vector::zero);

    Omega_ =   Axis_*omega_
               + mesh_.time().deltaT().value() * ( inv(I_) & T );

    //Split Omega_ into Axis_ and omega_
    omega_ = mag(Omega_);

    if(omega_ < 1e-32)
    {
        Axis_ = vector::zero;
    }
    else
    {
       Axis_ =  Omega_/omega_;
    }
}

//---------------------------------------------------------------------------//
void Foam::immersedBody::ensureStorageForCurrentMesh_()
{

    // Resize ONLY when size changed. Do NOT clear per-entry data here.


}

//------------------------------------------------------------------------//
//Create interpolation points
void
immersedBody::calculateInterpolationPoints
(
    volScalarField& body,
    triSurfaceSearch& ibTriSurfSearch
)
{
    double sqrtThree_ = sqrt(3.0);
    meshSearch search_(mesh_);

    //clear previous
    interpolationPoints_.clear();
    interpolationCells_.clear();

    //Create temporary surface normals
    volVectorField surfNorm(-fvc::grad(body));
    // If nothing on the surface yet, bail out cleanly
    // if (surfCells_.empty())
    // {
    //     WarningInFunction
    //         << "No surface cells at time " << mesh_.time().timeName()
    //         << ". Skipping interpolation-point build." << nl;
    //     interpolationCells_.clear();
    //     interpolationPoints_.clear();
    //     return;   // <-- important // NO, YOU CAN"T DO THIS
    // }

    //Create local list of remote nodes
    List<point> localRemotePoints;

    forAll(surfCells_,cell)
    {
        //get surface cell label
        label scell = surfCells_[cell];

        //create vector for points and cells and add to main vectors
        pointField intPoints;
        interpolationPoints_.append(intPoints);

        labelList intCells;
        interpolationCells_.append(intCells);

        //Get interpolation distance
	const scalar intDist = sqrtThree_*std::pow(mesh_.V()[scell], 1.0/3.0);

	// Build a safe outward direction (avoid divide-by-zero)
	vector intVec(vector::zero);
	const scalar nmag = mag(surfNorm[scell]);
	if (nmag > VSMALL)
	  {
	    intVec = (surfNorm[scell]/nmag) * intDist;
	  }
	else
	  {
	    // final fallback if the gradient is (near) zero
	    intVec = vector(0,0,1) * intDist;
	  }

        //Approximate distance using body
        point surfPoint = mesh_.C()[scell] + intVec*(0.5-body[scell]);

        //Add to list
        interpolationPoints_[cell].append(surfPoint);

       //Add other interpolation points
        for(int order=0;order<ORDER;order++)
        {
            surfPoint = surfPoint + intVec;
            interpolationPoints_[cell].append(surfPoint);
        }

        //Get cells
        for(int order=0;order<ORDER;order++)
        {

            label cellI;

            //Check if inside the domain (if not, then set to -1)
            if
            (
                !search_.isInside
                (
                    interpolationPoints_[cell][order+1]
                )
            )
            {
                cellI = -1;
                localRemotePoints.append(interpolationPoints_[cell][order+1]);
            }
            else
            {
                cellI =
                search_.findCell(interpolationPoints_[cell][order+1]);
            }

            interpolationCells_[cell].append(cellI);
        }

    }

    //- Parallel communication

    if (!Pstream::parRun())
    {
        return;
    }

    //- Standard MPI algorithm in OF language, nothing special.
    //  Uninterested reader can skip the code.
    labelList displ(UPstream::nProcs(),0);
    labelList numfrags(UPstream::nProcs(),0);
    label   remoteSize = 0;

    numfrags[UPstream::myProcNo()] = localRemotePoints.size();

    reduce(numfrags,sumOp<labelList>());

    forAll(displ,proc)
    {
        displ[proc] = remoteSize;
        remoteSize += numfrags[proc];
    }

    remoteDispl_ = displ[UPstream::myProcNo()];

    remotePoints_.clear();
    remotePoints_.resize(remoteSize,vector::zero);

    forAll(localRemotePoints,lrpI)
    {
        remotePoints_[remoteDispl_+lrpI] = localRemotePoints[lrpI];
    }

    reduce(remotePoints_,sumOp<List<point> >());

    //- Now every processor has a list of unfound (-1) intepolation points.
    //  At this point the algorithm establishes who holds what.

    //- Only poistion in remotePoints_ and corresponding cell are stored.
    remoteCells_.clear();

    forAll(remotePoints_,rpI)
    {
        //- Continue loop if outside the domain.
        //  Perhaps someone else has it.
        if(!search_.isInside(remotePoints_[rpI]))
        {
            continue;
        }

        //-  It should be here, run findCell
        label cellI =
        search_.findCell(remotePoints_[rpI]);

        //- Check if valid (it should be), then add to list
        if(cellI>-1)
        {
            labelList tmpList(2,zero());
            tmpList[0] = rpI;
            tmpList[1] = cellI;

            remoteCells_.append(tmpList);
        }

    }

}
//---------------------------------------------------------------------------//
void immersedBody::calculateGeometricalProperties( volScalarField& body )
{

     ensureStorageForCurrentMesh_();

    //Get density
    //dimensionedScalar rho(immersedDict_.lookup("rho"));
    const dimensionedScalar rho("rho", immersedDict_);

    //Evaluate center of mass
    M_ = 0.;
    vector tmpCom(vector::zero);
    // Preserve previous values in case M_==0 this step
    const vector prevCoM(CoM_);
    const symmTensor prevI(I_);
    // (kept) CoM_ left unchanged here; will update after mass check
    I_ = symmTensor::zero;

    // If no support cells at all, keep previous CoM/I and leave M_=0
    if (intCells_.empty() && surfCells_.empty())
    {
        WarningInFunction << "No body-support cells at time " << mesh_.time().timeName() << ". Keeping previous CoM/I." << nl;
        M_ = 0;

    }

    for(int cell=0;cell<intCells_.size()+surfCells_.size();cell++)
    {
        label cellI;

        if(cell<intCells_.size())
        {
            cellI = intCells_[cell];
        // Defend against any stale/out-of-range indices

        if (cellI < 0 || cellI >= mesh_.nCells())

        {

            WarningInFunction << "Bad cell index " << cellI

                << " (nCells=" << mesh_.nCells() << "), skipping." << nl;

            continue;

        }

        }
        else
        {
            cellI = surfCells_[cell-intCells_.size()];
        // Defend against any stale/out-of-range indices
        if (cellI < 0 || cellI >= mesh_.nCells())
        {
            WarningInFunction << "Bad cell index " << cellI
                << " (nCells=" << mesh_.nCells() << "), skipping." << nl;
            continue;
        }

        }

        M_ += body[cellI] * rho.value() * mesh_.V()[cellI];

        tmpCom  +=   body[cellI] * rho.value()
                   * mesh_.V()[cellI] * mesh_.C()[cellI];

        I_.xx() += body[cellI]*rho.value()*mesh_.V()[cellI]
                  * (
                        mesh_.C()[cellI].y()*mesh_.C()[cellI].y()
                      + mesh_.C()[cellI].z()*mesh_.C()[cellI].z()
                    );

        I_.yy() += body[cellI]*rho.value()*mesh_.V()[cellI]
                  * (
                       mesh_.C()[cellI].x()*mesh_.C()[cellI].x()
                     + mesh_.C()[cellI].z()*mesh_.C()[cellI].z()
                    );

       I_.zz() += body[cellI]*rho.value()*mesh_.V()[cellI]
                 * (
                        mesh_.C()[cellI].y()*mesh_.C()[cellI].y()
                      + mesh_.C()[cellI].x()*mesh_.C()[cellI].x()
                   );

        I_.xy() -= body[cellI]*rho.value()*mesh_.V()[cellI]
                 * (
                        mesh_.C()[cellI].x()*mesh_.C()[cellI].y()
                   );

        I_.xz() -= body[cellI]*rho.value()*mesh_.V()[cellI]
                 * (
                        mesh_.C()[cellI].x()*mesh_.C()[cellI].z()
                   );

        I_.yz() -= body[cellI]*rho.value()*mesh_.V()[cellI]
                 * (
                        mesh_.C()[cellI].y()*mesh_.C()[cellI].z()
                   );

    }

    //Collect from processors
    reduce(M_, sumOp<scalar>());
    reduce(tmpCom,  sumOp<vector>());
    reduce(I_,  sumOp<symmTensor>());

    if (mag(M_) > VSMALL) { CoM_ = tmpCom / M_; } else { WarningInFunction << "Zero immersed mass at time " << mesh_.time().timeName() << "; keeping previous CoM/I." << nl; CoM_ = prevCoM; I_ = prevI; }
}
//---------------------------------------------------------------------------//
//Move immersed body according to body operation
void immersedBody::moveImmersedBody()
{
    //Rotation angle
    scalar angle = omega_*mesh_.time().deltaT().value();
    vector transIncr = Vel_*mesh_.time().deltaT().value();

    // pointField bodyPoints (bodySurfMesh_->points());

    tmp<pointField> tBody = bodySurfMesh_->points();
    pointField bodyPoints(tBody());

    //- Just take one point to compute the reference direction
    vector n1(bodyPoints[0]) ;

    //- Create vector normal to axis and radius
    vector normV = (n1-CoM_)^Axis_;

    //Move in tangential direction
    point n2 = bodyPoints[0] +angle*normV;

    //- Rescale
    n1 /= mag(n1);
    n2 /= mag(n2);

    //- Rotate points from axis n1 to axis n2
    tensor T(rotationTensor(n1, n2));

    //- Rotate points
    bodyPoints = transform(T, bodyPoints);

    if (scaleEnabled_)
      {
	scalar t = mesh_.time().value();

	scalar s = s0_;
	if      (t <= t0_) s = s0_;
	else if (t >= t1_) s = s1_;
	else               s = s0_ + (s1_ - s0_)*(t - t0_)/max(SMALL, t1_ - t0_);

	// incremental factor
	const scalar sInc = s / max(VSMALL, sPrev_);
	forAll(bodyPoints, p)
	  {
	    bodyPoints[p] = CoM_ + sInc*(bodyPoints[p] - CoM_);
	  }
	sPrev_ = s;
      }

    //- Translate points
    forAll(bodyPoints,p)
    {
        bodyPoints[p] += transIncr;
    }

    //move mesh
    bodySurfMesh_->movePoints(bodyPoints);

}
//---------------------------------------------------------------------------//
//Reset body field for this immersed object
void immersedBody::resetBody(volScalarField& body)
{
    ensureStorageForCurrentMesh_();
    // Mesh may have changed; drop any stale selections and reset the indicator field
    body = dimensionedScalar("zero", dimless, 0.0);
    intCells_.clear();
    surfCells_.clear();
    surfNorm_.clear();

    //Simply loop over all the cells and set to zero
    for(int cell=0;cell<intCells_.size()+surfCells_.size();cell++)
    {
        label cellI;

        if(cell<intCells_.size())
        {
            cellI = intCells_[cell];
        // Defend against any stale/out-of-range indices

        if (cellI < 0 || cellI >= mesh_.nCells())

        {

            WarningInFunction << "Bad cell index " << cellI

                << " (nCells=" << mesh_.nCells() << "), skipping." << nl;

            continue;

        }

        }
        else
        {
            cellI = surfCells_[cell-intCells_.size()];
        // Defend against any stale/out-of-range indices
        if (cellI < 0 || cellI >= mesh_.nCells())
        {
            WarningInFunction << "Bad cell index " << cellI
                << " (nCells=" << mesh_.nCells() << "), skipping." << nl;
            continue;
        }

        }

        body[cellI] = 0.0;
    }

}
//---------------------------------------------------------------------------//
//Refine body field for this immersed object using MC-like algorithm
//Cells are assumed to be hexahedral at the particle surface
//(but can have different edge length)
void immersedBody::refineBody
(
    volScalarField& body,
    triSurfaceSearch& ibTriSurfSearch,
    const pointField& pp
)
{
    if(!immersedDict_.found("refineMC"))
    {
        return;
    }


    scalar nPointsEdge = readScalar(immersedDict_.lookup("refineMC"));

    //loop over all the surface cells
    forAll(surfCells_,cell)
    {

        label cellI = surfCells_[cell];

        scalar deltaV = 1.0/(nPointsEdge*nPointsEdge*nPointsEdge);

        //Get cell center
        point centerC = mesh_.C()[cellI];

        //Get one node
        //Check if partially or completely inside
        const labelList& vertexLabels = mesh_.cellPoints()[cellI];
        const pointField vertexPoints(pp,vertexLabels);
        point baseNode = vertexPoints[0];

        //create vector representing 3d diagonal of the cell
        vector edgesC = 2.0*(centerC - baseNode);

        //create list of points
        pointField pointsMC;

        //create deltas
        scalar delta_i = edgesC[0]/nPointsEdge;
        scalar delta_j = edgesC[1]/nPointsEdge;
        scalar delta_k = edgesC[2]/nPointsEdge;

        //add points to list
        for(int i=0;i<nPointsEdge;i++)
        {
            //point i-coordinate
            scalar icoord = baseNode[0] + delta_i*(i+0.5);

            for(int j=0;j<nPointsEdge;j++)
            {
                //point j-coordinate
                scalar jcoord = baseNode[1] + delta_j*(j+0.5);

                for(int k=0;k<nPointsEdge;k++)
                {
                    //point k-coordinate
                    scalar kcoord = baseNode[2] + delta_k*(k+0.5);

                    //create point
                    point p(icoord,jcoord,kcoord);

                    //add to list
                    pointsMC.append(p);

                }
            }
        }

        //Check who is inside
        pointField pField(pointsMC);
        boolList pInside = ibTriSurfSearch.calcInside( pField );

        //Calculate new body
        scalar newbody = 0.0;

        forAll(pInside,p)
        {
            if(pInside[p])
            {
                newbody+=deltaV;
            }

        }


        body[cellI] = newbody;

    }

}
//---------------------------------------------------------------------------//
