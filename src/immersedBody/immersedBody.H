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

  Description
  class for immersed bodies representation.

  SourceFiles
  immersedBodies.C

  Contributors
  Federico Municchi (2018)
  \*---------------------------------------------------------------------------*/

#ifndef immersedBody_H
#define immersedBody_H

#include "dictionary.H"
#include "coordinateSystem.H"
#include "coordinateSystems.H"
#include "wordList.H"
#include "labelList.H"
#include "dimensionedScalar.H"
#include "dimensionedTensor.H"
#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "fvMatricesFwd.H"
#include "fvMesh.H"
#include "triSurfaceMesh.H"
#include "polyMesh.H"
#include "fvCFD.H"
#include "List.H"
#include "autoPtr.H"
#include "UPstream.H"
#include "interpolationCellPoint.H"
#include "interpolationCell.H"
#include "triSurfaceSearch.H"


namespace Foam
{
  enum {
    STATICBODY, //Body does not move
    TRANSROTATINGBODY, //Trans-rotational imposed motion
    FLUIDCOUPLING  //Fluid driven motion
  };

  class immersedBody
  {

  private:

    //- Name of immersed body
    word          bodyName_;

    //- First time update?
    bool          isFirstUpdate_;

    //- Dictionary
    dictionary&   immersedDict_;

    //- Mesh
    const fvMesh& mesh_;

    //- Dictionary
    dictionary&   transportProperties_;


    //- Body mass
    scalar M_;

    //- Location of center of mass
    vector CoM_;

    //- Rotation axis
    vector Axis_;

    //- Angular velocity
    scalar omega_;

    //- Body velocity
    vector Vel_;

    // --- deformation (optional scaling over time) ---
    bool   scaleEnabled_;      // on/off
    scalar t0_;                // start time of scaling
    scalar t1_;                // end time of scaling
    scalar s0_;                // scale at t0 (1 = original size)
    scalar s1_;                // scale at t1
    // private:
    scalar sPrev_;   // remember last absolute scale


    // ===== leafletRT: prescribed z(r,t) motion & autoR0 (added) =====

    // ===== leafletRT: prescribed z(r,t) motion & autoR0 (added) =====

    bool   leafletEnabled_ = false;  // on when deformingBody{type leafletRT;} is parsed
    vector leafletAxis_;             // unit axis for axial motion (from dict "axis")
    point  leafletCenter_;           // origin for r (from dict "center")

    // Rim radius and clamp (default-initialised; validated in .C)
    scalar R0_           = 0.0;      // outer radius; computed if autoR0_ = true
    scalar Rmin_         = 0.0;      // inner clamp radius (usually 0)
    scalar Rmax_         = 0.0;      // outer clamp radius (0 => use R0_ at runtime)
    bool   clampOutside_ = true;     // freeze points with R outside [Rmin_, Rcap]
    bool   autoR0_       = false;    // compute R0_ from STL points

    // Tables for radial law f(r̂), amplitude A(t), and baseline B(t)
    List<scalar> rHat_, fR_;         // r̂ in [0..1], f(r̂)
    List<scalar> tA_,   A_;          // time vs amplitude
    List<scalar> tB_,   B_;          // time vs offset

    // ===== end leafletRT additions ====
    // ===== leafletPlane: prescribed Δ = A(t)*g(d̂) + B(t) along 'axis', d̂ ∈ [0..1] =====

    bool   planeEnabled_ = false;   // on when deformingBody{type leafletPlane;} is parsed
    vector planeNormal_;            // (unit) normal of the hinge plane
    point  planePoint_;             // a point on the hinge plane (e.g. (0,0.05,0))
    bool   autoH0_      = false;    // compute H0_ from STL extents along planeNormal
    scalar H0_          = 0.0;      // span from hinge plane to tip (metres)

    //========================================//

    //- Moment of intertia tensor
    symmTensor  I_;

    //- Operation to perform
    int   bodyOperation_;

    //- Base geometry
    autoPtr<triSurfaceMesh>             bodySurfMesh_;

    //- Surface cells
    labelList           surfCells_;

    //- Internal cells
    labelList            intCells_;

    //- Surface normals
    vectorField         surfNorm_;

    //- Interpolation points, point[0] is a surface point
    List<pointField>     interpolationPoints_;

    //- Interpolation cells
    List<labelList>              interpolationCells_;

    //- List of remote points (global)
    List<point>        remotePoints_;

    //- List of remote cells
    List<labelList>    remoteCells_;

    //- Processor displacement in remotePoints_
    scalar      remoteDispl_;


    //- Create immersed body info
    void  createImmersedBody( volScalarField& body );

    //- Update immersed body info
    void  updateImmersedBody
    (
     volScalarField& body,
     volVectorField& f
     );

    //- Create interpolation points
    void  calculateInterpolationPoints
    (
     volScalarField& body,
     triSurfaceSearch& ibTriSurfSearch
     );

    //- Evaluate center of mass and inertia tensor
    void calculateGeometricalProperties(volScalarField& body);

    //- Rotate immersed body
    void moveImmersedBody();

    //- Reset body field
    void resetBody(volScalarField& body);

    // Ensure internal storage matches current mesh sizes (after refine/unrefine)
    void ensureStorageForCurrentMesh_();

    //- Refine body field using MC-like integration
    void refineBody
    (
     volScalarField& body,
     triSurfaceSearch& ibTriSurfSearch,
     const pointField& pp
     );

    //- Update Fluid-Solid coupling
    void updateCoupling
    (
     volScalarField & body,
     volVectorField & f
      );

    //- Tranform body using dictionary
    void transformBody(dictionary& transformDict);

    // ===== leafletRT helpers (declared; definitions in .C) =====

    // piecewise-linear interpolation helper
    scalar interp(const List<scalar>& x, const List<scalar>& y, scalar xv) const;

    // sampled amplitude/offset in time
    inline scalar AofT(const scalar t) const { return interp(tA_, A_, t); }
    inline scalar BofT(const scalar t) const { return interp(tB_, B_, t); }

    // sampled radial shape f(R/R0)
    inline scalar fOfR(const scalar R) const
    {
      const scalar rhat = (R0_ > SMALL ? Foam::max(0.0, Foam::min(1.0, R/R0_)) : 0.0);
      return interp(rHat_, fR_, rhat);
    }

    // compute R0_ as the max radius from STL points in plane ⟂ leafletAxis_
    void computeAutoR0_();

    // read "deformingBody { type leafletRT; ... }" sub-dictionary
    void readLeafletRT(const dictionary& d);

    // apply prescribed motion to geometry points for current time
    void applyLeafletRT();

    // leafletPlane helpers
    void computeAutoH0_();                       // find max |distance| to the plane
    void readLeafletPlane(const dictionary& d);  // parse deformingBody{…}
    void applyLeafletPlane();                    // move points for current time

    // ===== end leafletRT helpers ====

    // ===== annulusRZ: prescribed ri(z,t)/ro(z,t) motion in cylindrical coords =====
    bool   annulusEnabled_ = false;   // on when deformingBody{type annulusRZ;} is parsed
    vector annAxis_;                  // unit z-axis of leaflet (from dict "axis")
    point  annCenter_;                // origin for r,z (from dict "center")

    scalar ri0_   = 0.0;              // initial inner radius
    scalar zf_    = 0.0;              // flatten (saturation) height
    scalar thick_ = 1.0;              // thickness (ro = ri + thick)
    scalar T_     = 1.0;              // period for s(t) = sin^2(pi t / T)
    scalar tol_   = 1e-6;             // band classification tolerance (radial)
    scalar wallTol_ = 4e-4;           // ± window around each wall that will move
    scalar p_       = 3.0;            // concavity exponent for ease-out in z (p=3 recommended)

    // NEW: base-clamp & rest-geometry bookkeeping
    scalar z0_ = 0.0;                 // base plane z (rest STL, along annAxis_)
    scalar L_  = 0.0;                 // total length (zMax - zMin) from rest STL
    scalar tolZ_ = 1e-6;              // axial tolerance for base/top ring tagging
    word   bandPolicy_{"baseFixed"};  // "baseFixed" (clamp base) | "fixed" (legacy)

    // Per-vertex classification on the REST STL
    //   0 = inner skin, 1 = band, 2 = outer skin, -1 = untouched
    List<label>  pointClass_;
    // For band points: α in [0,1] = (R0 - ri0)/thickness from the REST geometry
    List<scalar> alpha_;
    // Base/top rings on REST STL
    List<bool>   isBase_, isTop_;
    // Snapshot of REST points to avoid cumulative drift
    pointField   X0_;

    // Helpers (defined in .C)
    scalar sOfT(const scalar t) const;                     // sin^2(pi t / T)
    scalar uOfZ(const scalar z) const;                     // clamp((z - z0_)/zf,0,1)
    scalar riOfZT(const scalar z, const scalar t) const;   // max(ri0*(1 - s*u^2), 0)
    scalar roOfZT(const scalar z, const scalar t) const;   // ri + thick

							      // API
    void readAnnulusRZ(const dictionary& d);
    void applyAnnulusRZ();

    // ===== Mitral leaflet: inner-surface axis + per-slice radial closure =====
    bool   mitralEnabled_ = false;
    // Optional: flip "top" cap normal toward a reference valve STL
    Switch  m_flipNormalTowardValve_{false};
    fileName m_valveRefFile_ = fileName("");

    // Geometry (inner-surface based)
    point  mCt_;                 // top-cap (annulus) center on INNER surface
    point  mCb_;                 // bottom-cap (free edge) center on INNER surface
    vector mAxisUnit_;           // unit axis from Ct -> Cb
    scalar mLen_;                // L = |Cb - Ct|

    // Time & space law
    scalar mPeriod_    = 1.0;
    scalar mNearRamp_  = 0.004;
    scalar mActiveLen_ = 0.016;
    scalar mFmax_      = 0.995;
    word   mTimeLaw_   = word("cos2");
    word   mSpaceLaw_  = word("smootherstep");   // NEW: C2 by default

    // REST snapshot (inner surface points only)
    pointField mX0_;              // original points (from triSurfaceMesh)
    // Per-vertex bookkeeping
    List<scalar> mXi_;            // axial coordinate ξ in [0,L]
    List<vector> mRperp0_;        // in-plane vector at t=0 (X0 - C(ξ)) projected ⟂ axis
    List<scalar> mR0_;            // radius r = |r_perp0|
    List<scalar> mAlpha_;         // through-thickness fraction [0..1] if paired; else 0
    List<vector> mNthru_;         // through-thickness unit vector (if available), else (0,0,0)

    // If you have paired inner/outer meshes, set tLeaf_ > 0; else keep 0 and no thickness enforcement.
    scalar tLeaf_ = 0.0;          // leaflet thickness [m]; 0 => single-skin

    // Reading & application
    void readMitralSliceAxis(const dictionary& d);
    void applyMitralSliceAxis();

    // Utilities (shared)
    void autoFitPlaneFromSurface(vector& nHat, point& centroid) const;   // (you already have this)
    void maybeFlipNormalUsingValve(vector& nHat, const point& centroid) const; // (you already have this)
    void fitTopBottomCapsFromSurface(point& Ct, point& Cb, vector& axisUnit) const;

    // Time gains
    inline scalar timeGain_(scalar t) const
    {
        scalar x; // normalized 0..1 in period
        if (mPeriod_ > VSMALL)
        {
            scalar tCycle = std::fmod(t, mPeriod_);
            if (tCycle < 0) tCycle += mPeriod_;
            x = tCycle / mPeriod_;
        }
        else x = 0;

        if (mTimeLaw_ == "cos2")        return 0.5*(1.0 - Foam::cos(2.0*Foam::constant::mathematical::pi*x));
        else if (mTimeLaw_ == "smoothstep")
        {
            scalar s = clamp(x, 0.0, 1.0);
            return (3*s*s - 2*s*s*s);
        }
        // default sin^2
        const scalar s = Foam::sin(Foam::constant::mathematical::pi * x);
        return s*s;
    }

    inline scalar spaceGain_(scalar xi) const
    {
      if (xi < 0) return 0.0;
      if (xi <= mNearRamp_)
	{
	  const scalar s = clamp(xi / max(VSMALL, mNearRamp_), 0.0, 1.0);
	  if (mSpaceLaw_ == "linear")
            return s;
	  else if (mSpaceLaw_ == "smoothstep")        // C1
            return (3*s*s - 2*s*s*s);
	  else /* "smootherstep" (default) */          // C2
            return (s*s*s * (10.0 - 15.0*s + 6.0*s*s));
	}
      if (xi <= mActiveLen_) return 1.0;    // derivative 0 from left; constant right
      return 0.0;
    }


  public:

    immersedBody
    (
     word fileName,
     const Foam::fvMesh& mesh,
     dictionary& HFDIBDict,
     dictionary& transportProperties
     );

    ~immersedBody();

    const labelList& getSurfaceCellList() const
    {
      return surfCells_;
    };

    const List<pointField>& getInterpolationPoints() const
    {
      return interpolationPoints_;
    };

    const List<labelList>& getInterpolationCells() const
    {
      return  interpolationCells_;
    };

    //- Update body field
    void updateBodyField
    (
     volScalarField & body,
     volVectorField & f
     );

    //- Update imposed field
    template<class Type>
    void updateImposedField(Type&, word);

    //- Remote interpolation
    template<class Type>
    void remoteInterpolateField
    (
     autoPtr<interpolation<Type> >&  interpField,
     List<Type>&  interpValues
     );

  };

#   include "immersedBodyTemplates.C"

}
#endif
