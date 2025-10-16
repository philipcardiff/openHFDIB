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
#include "pointIndexHit.H"
#include "mathematicalConstants.H"
#include <cctype>   // for std::tolower



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
I_(symmTensor::zero),
bodySurfMesh_
(
    new triSurfaceMesh
    (
        IOobject
        (
            fileName(immersedDict_.lookup("fileName")),
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
    Info<< endl;

    if(immersedDict_.found("transform"))
    {
        dictionary transformDict = immersedDict_.subDict("transform");

        transformBody(transformDict);
    }

    //Return if declared static
    if(immersedDict_.found("staticBody") )
    {
        bodyOperation_=STATICBODY;
        Info<< bodyName << " is static body." << endl;
    }
    else if(immersedDict_.found("transRotatingBody"))
    {
        bodyOperation_=TRANSROTATINGBODY;

        //Get basic quantities from dict
        Axis_ =
            vector(immersedDict_.subDict("transRotatingBody").lookup("axis"));
        CoM_ =
            vector(immersedDict_.subDict("transRotatingBody").lookup("center"));

        omega_  = readScalar
        (
            immersedDict_.subDict("transRotatingBody").lookup("omega")
        );

        Vel_ =
            vector(immersedDict_.subDict("transRotatingBody").lookup("velocity"));

        Info<< bodyName << " has scripted trans rotational motion." << endl;
    }
    // --- enable mitral slice-axis deformation if present ---
    else if (immersedDict_.found("deformingBody_mitral"))
    {
        const dictionary& d = immersedDict_.subDict("deformingBody_mitral");
        readMitralSliceAxis(d);
        mitralEnabled_ = true;
        Info<< bodyName_ << " : deformingBody_mitral (slice-axis) enabled" << nl;
    }
    else if (immersedDict_.found("fluidCoupling"))
    {
        bodyOperation_ = FLUIDCOUPLING;
        Info<< bodyName << " is coupled with fluid phase." << endl;

    }
    else
    {
        Info<< "No body operation was found for " << bodyName << endl
             << "Assuming static body.";
        bodyOperation_=STATICBODY;
    }
}
//---------------------------------------------------------------------------//
immersedBody::~immersedBody()
{
    bodySurfMesh_.clear();
}
//---------------------------------------------------------------------------//
void immersedBody::transformBody(dictionary& transformDict)
{

    Info<< "Transforming immersed body " << bodyName_
         << " using dictionary" << endl;

    pointField bodyPoints(bodySurfMesh_->points());

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
void immersedBody::updateBodyField
(
    volScalarField& body,
    volVectorField& f
)
{
    if (isFirstUpdate_)
    {
        createImmersedBody(body);
        isFirstUpdate_ = false;
    }
    else
    {
        updateImmersedBody(body, f);
    }

    body.correctBoundaryConditions();
}

//---------------------------------------------------------------------------//
//Create immersed body info
void immersedBody::createImmersedBody(volScalarField& body)
{

    triSurface ibTemp(bodySurfMesh_());
    triSurfaceSearch ibTriSurfSearch(ibTemp);
    const pointField & pp = mesh_.points();

    intCells_.clear();
    surfCells_.clear();
    surfNorm_.clear();

    //Info<<"\n Creating immersed body\n";

    //Fill body field with first estimation
    forAll(mesh_.C(),cellI)
    {
        //Check if partially or completely inside
        const labelList& vertexLabels = mesh_.cellPoints()[cellI];
        const pointField vertexPoints(pp,vertexLabels);
        boolList vertexesInside = ibTriSurfSearch.calcInside(vertexPoints);

        bool bodyCell(false);
        body[cellI] = 0.;

        vector centerOut(vector::zero);

        forAll(vertexesInside, verIn)
        {
            if(vertexesInside[verIn]==true)
            {
                //fraction of cell covered
                body[cellI] += 1.0/(vertexPoints.size());
               // Info<< "Found vertex inside\n";
                bodyCell = true;
             }

        }

        //Add to corresponding vector
        if(bodyCell)
        {
            if( body[cellI]>0.9)
            {
                intCells_.append(cellI);
            }
            else if (body[cellI]>0.1)
            {
                surfCells_.append(cellI);

            }
        }
    }

    //refine body as stated in the dictionary
    refineBody(body,ibTriSurfSearch,pp);
    calculateInterpolationPoints(body,ibTriSurfSearch);
    calculateGeometricalProperties(body);

}
//---------------------------------------------------------------------------//
//Update immersed body info
void immersedBody::updateImmersedBody
(
    volScalarField & body,
    volVectorField & f
)
{
    //Check Operation to perform

    if (bodyOperation_==STATICBODY && !mitralEnabled_)
    {
        return;
    }
    else
    {
        if (bodyOperation_==FLUIDCOUPLING)
        {
            updateCoupling(body,f);
        }

        // Geometry-dependent fields must be recomputed after motion
        resetBody(body);

        // Apply mitral deformation if enabled
        if (mitralEnabled_)
	{
	  applyMitralSliceAxis();
	}

        // Legacy trans-rotational motion if requested
        if (bodyOperation_==TRANSROTATINGBODY)
	{
	    moveImmersedBody();
	 }
    }

    //TODO:inefficient algorithm...
    createImmersedBody(body);

}
//---------------------------------------------------------------------------//
void immersedBody::updateCoupling
(
    volScalarField & body,
    volVectorField & f
)
{

    const uniformDimensionedVectorField& g =
        mesh_.lookupObject<uniformDimensionedVectorField>("g");

    const dimensionedScalar rhof("rho", transportProperties_);

    const dimensionedScalar rho("rho", immersedDict_);

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
                F / M_
                + (1.0-rhof.value()/rho.value())*g.value()
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
//----------------------------------------------------------------------------//
void immersedBody::readMitralSliceAxis(const dictionary& d)
{
    // Period & gains
    mPeriod_    = d.lookupOrDefault<scalar>("period",    1.0);
    mNearRamp_  = d.lookupOrDefault<scalar>("nearRamp",  0.004);
    mActiveLen_ = d.lookupOrDefault<scalar>("activeLen", 0.016);
    mFmax_      = d.lookupOrDefault<scalar>("Fmax",      0.995);
    mTimeLaw_   = d.lookupOrDefault<word>  ("timeLaw",   word("cos2"));
    tLeaf_      = d.lookupOrDefault<scalar>("leafThickness", 0.0);
    mSpaceLaw_  = d.lookupOrDefault<word>("spaceLaw", word("smootherstep"));

     // NEW: which end is fixed for spatial ramp (default: "top")
    {
        word fixedEndWord = d.lookupOrDefault<word>("fixedEnd", word("top"));
	// fixedEndWord = fixedEndWord.tolower();
	// portable lowercase (Foam::word has no .tolower())
	Foam::string fe = fixedEndWord;          // copy to a mutable string
	for (char& ch : fe) {
	  ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
	}
	fixedEndWord = Foam::word(fe);

        if (fixedEndWord == "bottom" || fixedEndWord == "cb")
        {
            mFixedEnd_ = BottomFixed;
        }
        else
        {
            // accept "top", "ct", or anything else -> TopFixed
            mFixedEnd_ = TopFixed;
        }
    }

    // Optional orientation helpers
    m_flipNormalTowardValve_ = d.lookupOrDefault<Switch>("flipNormalTowardValve", false);
    m_valveRefFile_          = d.lookupOrDefault<fileName>("valveRefFile", fileName(""));

    // Optional timing windows
    if (mTimeLaw_ == "twoWindow")
    {
        if (d.found("closeWindow")) { const List<scalar> cw(d.lookup("closeWindow")); if (cw.size()==2){ mClose0_=cw[0]; mClose1_=cw[1]; } }
        if (d.found("openWindow"))  { const List<scalar> ow(d.lookup("openWindow"));  if (ow.size()==2){ mOpen0_=ow[0];  mOpen1_=ow[1];  } }
    }
    else if (mTimeLaw_ == "threeWindow")
    {
        if (d.found("close1Window")) { const List<scalar> c1(d.lookup("close1Window")); if (c1.size()==2){ mClose0_=c1[0]; mClose1_=c1[1]; } }
        if (d.found("openWindow"))   { const List<scalar> ow(d.lookup("openWindow"));   if (ow.size()==2){ mOpen0_=ow[0];  mOpen1_=ow[1];  } }
        if (d.found("close2Window")) { const List<scalar> c2(d.lookup("close2Window")); if (c2.size()==2){ mClose2_0_=c2[0]; mClose2_1_=c2[1]; } }
    }

    // Clamp all window endpoints to [0,1]
    auto clamp01 = [](scalar v)->scalar { return max(0.0, min(1.0, v)); };
    mClose0_   = clamp01(mClose0_);   mClose1_   = clamp01(mClose1_);
    mOpen0_    = clamp01(mOpen0_);    mOpen1_    = clamp01(mOpen1_);
    mClose2_0_ = clamp01(mClose2_0_); mClose2_1_ = clamp01(mClose2_1_);

    // Optional explicit Ct/Cb; otherwise auto-detect from STL
    const bool haveCt = d.found("topCenter");
    const bool haveCb = d.found("bottomCenter");
    if (haveCt) mCt_ = point(d.lookup("topCenter"));
    if (haveCb) mCb_ = point(d.lookup("bottomCenter"));

    // Snapshot REST points from current triSurfaceMesh
    tmp<pointField> tPts = bodySurfMesh_->points();
    mX0_ = tPts();

    // Axis detection
    if (!haveCt || !haveCb) { fitTopBottomCapsFromSurface(mCt_, mCb_, mAxisUnit_); }
    else
    {
        mAxisUnit_ = vector(mCb_ - mCt_);
        const scalar L = mag(mAxisUnit_);
        mAxisUnit_ = (L > VSMALL ? mAxisUnit_/L : vector(0,0,1));
    }

    mLen_ = mag(vector(mCb_ - mCt_));
    if (mLen_ <= VSMALL)
    {
        FatalErrorInFunction << "Mitral slice-axis: Ct and Cb invalid; cannot define axis." << exit(FatalError);
    }

    // Precompute per-vertex ξ and in-plane rest vectors
    const label n = mX0_.size();
    mXi_.setSize(n, 0.0);
    mRperp0_.setSize(n, vector::zero);
    mR0_.setSize(n, 0.0);

    forAll(mX0_, i)
    {
        const vector d0 = vector(mX0_[i] - mCt_);
        scalar xi = (d0 & mAxisUnit_);
        xi = clamp(xi, 0.0, mLen_);
        mXi_[i] = xi;

        const vector Cxi = vector(mCt_) + xi*mAxisUnit_;
        vector r = vector(mX0_[i] - Cxi);
        r -= (r & mAxisUnit_) * mAxisUnit_;
        mRperp0_[i] = r;
        mR0_[i] = mag(r);
    }
}

void immersedBody::applyMitralSliceAxis()
{
    if (!mitralEnabled_) return;

    const scalar tNow = mesh_.time().value();
    const scalar gT   = timeGain_(tNow);   // 0..1

    // Rebuild from REST each step to avoid drift
    pointField pts(mX0_);

    label moved = 0; scalar minF = GREAT, maxF = -GREAT;

    forAll(pts, i)
    {
        const scalar xi = mXi_[i];
	// const scalar gS = spaceGain_(xi);
        // NEW: bottom-fixed uses distance from bottom cap instead of top
        const scalar xiLocal = (mFixedEnd_ == BottomFixed ? (mLen_ - xi) : xi);
        const scalar gS = spaceGain_(xiLocal);

        const scalar F  = clamp(mFmax_ * gS * gT, 0.0, 0.999);

        const vector r0 = mRperp0_[i];
        if (mag(r0) > VSMALL && F > SMALL)
        {
            const vector Cxi = vector(mCt_) + xi*mAxisUnit_;
            const vector Xnew = Cxi + (1.0 - F) * r0;
            pts[i] = point(Xnew);
            ++moved; minF = min(minF, F); maxF = max(maxF, F);
        }
        // else: keep REST
    }

    if ((mesh_.time().timeIndex() % 10) == 0)
    {
        if (moved == 0) { minF = 0; maxF = 0; }
        Info<< "mitral(slice-axis): t=" << tNow
            << " timeGain=" << gT
            << " movedPts=" << moved
            << " F[min,max]=[" << minF << "," << maxF << "]" << nl;
    }

    bodySurfMesh_->movePoints(pts);
}

// --- Helpers to infer axis & orient it consistently ---

void immersedBody::autoFitPlaneFromSurface(vector& nHat, point& centroid) const
{
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

    // Pick the flattest, large-area cluster as “cap”
    label  best = -1; scalar bestScore = -GREAT;
    forAll(clusterN, k)
    {
        scalar cosMin = 1.0;
        for (label j : clusterFaces[k]) cosMin = min(cosMin, clusterN[k] & rec[j].n);
        const scalar score = clusterA[k] * max(0.0, cosMin - cosTol);
        if (score > bestScore) { bestScore = score; best = k; }
    }
    auto clusterCentroid = [&](label k)->point
    {
        point csum = point::zero; scalar At = 0.0;
        for (label j : clusterFaces[k]) { csum += rec[j].A * rec[j].c; At += rec[j].A; }
        return (At > VSMALL ? csum/At : rec[clusterFaces[k].first()].c);
    };
    centroid = clusterCentroid(best);
    const vector nChosen = clusterN[best];
    const scalar nmag = mag(nChosen);
    nHat = (nmag > VSMALL ? nChosen / nmag : vector(0,0,1));
}

void immersedBody::maybeFlipNormalUsingValve(vector& nHat, const point& centroid) const
{
    if (!m_flipNormalTowardValve_ || m_valveRefFile_.empty()) return;

    triSurface tsRef(m_valveRefFile_);
    if (!tsRef.size()) return;

    point cRef = point::zero; scalar A = 0.0;
    forAll(tsRef, fi)
    {
        const labelledTri& tri = tsRef[fi];
        const point& a = tsRef.points()[tri[0]];
        const point& b = tsRef.points()[tri[1]];
        const point& c = tsRef.points()[tri[2]];
        const vector n = triPointRef(a,b,c).areaNormal();
        const scalar dA = 0.5*mag(n);
        const point cTri( (a+b+c)/3.0 );
        cRef += dA * cTri; A += dA;
    }
    if (A > VSMALL) cRef /= A; else return;

    const vector v = vector(cRef - centroid);
    if ( (v & nHat) < 0 ) nHat = -nHat;
}

void immersedBody::fitTopBottomCapsFromSurface(point& Ct, point& Cb, vector& axisUnit) const
{
    const triSurface& ts = bodySurfMesh_->surface();
    tmp<pointField> tP = bodySurfMesh_->points();
    const pointField& P = tP();

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

    // Pick two near-opposite large-area clusters
    label bestTop=-1, bestBot=-1; scalar bestArea = -GREAT;
    forAll(clusterN, k)
    {
        const vector nk = (mag(clusterN[k])>VSMALL ? clusterN[k]/mag(clusterN[k]) : vector(0,0,1));
        forAll(clusterN, j)
        {
            if (j==k) continue;
            const vector nj = (mag(clusterN[j])>VSMALL ? clusterN[j]/mag(clusterN[j]) : vector(0,0,1));
            const scalar c = nk & nj;
            if (c < -0.95)  // ~>171°
            {
                const scalar comb = clusterA[k] + clusterA[j];
                if (comb > bestArea) { bestArea = comb; bestTop = k; bestBot = j; }
            }
        }
    }
    if (bestTop<0 || bestBot<0)
    {
        // fallback to two largest clusters
        List<label> idx(clusterA.size());
        forAll(idx,i) idx[i]=i;
        std::sort(idx.begin(), idx.end(), [&](label i, label j){return clusterA[i]>clusterA[j];});
        bestTop = (idx.size()>0 ? idx[0] : -1);
        bestBot = (idx.size()>1 ? idx[1] : -1);
    }
    if (bestTop<0 || bestBot<0)
    {
        WarningInFunction<< "fitTopBottomCapsFromSurface: could not find two cap clusters; using defaults." << nl;
        Ct = point::zero; Cb = point(0,0,1); axisUnit = vector(0,0,1); return;
    }

    Ct = clusterCentroid(bestTop);
    Cb = clusterCentroid(bestBot);

    // Orient consistently using (optional) valve reference
    vector nTop = (mag(clusterN[bestTop])>VSMALL ? clusterN[bestTop]/mag(clusterN[bestTop]) : vector(0,0,1));
    maybeFlipNormalUsingValve(nTop, Ct);
    if ( (vector(Cb - Ct) & nTop) < 0 ) { point tmp=Ct; Ct=Cb; Cb=tmp; }

    axisUnit = vector(Cb - Ct);
    const scalar L = mag(axisUnit);
    axisUnit = (L>VSMALL ? axisUnit/L : vector(0,0,1));
}



//---------------------------------------------------------------------------//
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
        scalar intDist = sqrtThree_*std::pow( mesh_.V()[scell] , 1.0/3.0 );
         vector intVec = surfNorm[scell] * (intDist)/(mag(surfNorm[scell]));

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

    if (!UPstream::parRun())
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


void immersedBody::calculateGeometricalProperties(volScalarField& body)
{
    // Get density
    const dimensionedScalar rho("rho", immersedDict_);

    // Evaluate center of mass
    M_ = 0.;
    vector tmpCom(vector::zero);
    CoM_ = vector::zero;
    I_ = symmTensor::zero;

    for(int cell=0;cell<intCells_.size()+surfCells_.size();cell++)
    {
        label cellI;

        if(cell<intCells_.size())
        {
            cellI = intCells_[cell];
        }
        else
        {
            cellI = surfCells_[cell-intCells_.size()];
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

    CoM_ = tmpCom / M_;
}
//---------------------------------------------------------------------------//
//Move immersed body according to body operation
void immersedBody::moveImmersedBody()
{


    //Rotation angle
    scalar angle = omega_*mesh_.time().deltaT().value();
    vector transIncr = Vel_*mesh_.time().deltaT().value();

    pointField bodyPoints (bodySurfMesh_->points());

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

    //Simply loop over all the cells and set to zero
    for(int cell=0;cell<intCells_.size()+surfCells_.size();cell++)
    {
        label cellI;

        if(cell<intCells_.size())
        {
            cellI = intCells_[cell];
        }
        else
        {
            cellI = surfCells_[cell-intCells_.size()];
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
