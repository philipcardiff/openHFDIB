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
#include "labelList.H"
#include "UPstream.H"
#include "interpolationCellPoint.H"
#include "interpolationCell.H"
#include <cmath>
#include "mathematicalConstants.H"

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

         // ===== Mitral leaflet: inner-surface axis + per-slice radial closure =====
        bool   mitralEnabled_ = false;
        // Optional orientation helper (can be left default)
        Switch  m_flipNormalTowardValve_{false};
        fileName m_valveRefFile_ = fileName("");

        // Geometry (inner-surface based)
        point  mCt_;                 // top-cap (annulus) center on INNER surface
        point  mCb_;                 // bottom-cap (free edge) center on INNER surface
        vector mAxisUnit_;           // unit axis from Ct -> Cb
        scalar mLen_;                // L = |Cb - Ct|

        // Time & space laws
        scalar mPeriod_    = 1.0;
        scalar mNearRamp_  = 0.004;
        scalar mActiveLen_ = 0.016;
        scalar mFmax_      = 0.995;
        word   mTimeLaw_   = word("cos2");
        scalar mClose0_{0.625}, mClose1_{0.655};
        scalar mOpen0_{0.940},  mOpen1_{1.000};
        scalar mClose2_0_{0.0}, mClose2_1_{0.0};   // optional second close window
        word   mSpaceLaw_  = word("smootherstep");

       // REST snapshot & per-vertex bookkeeping
       pointField   mX0_;           // rest points
       List<scalar> mXi_;           // axial coordinate ξ∈[0,L]
       List<vector> mRperp0_;       // in-plane vector at t=0 (X0 - C(ξ)) projected ⟂ axis
       List<scalar> mR0_;           // |mRperp0_|
       scalar       tLeaf_ = 0.0;   // leaflet thickness [m]; 0 => single-skin
       enum FixedEnd { TopFixed, BottomFixed };
       FixedEnd mFixedEnd_ = TopFixed;  // default matches prior mitral behavior

          // Reading & application
    void readMitralSliceAxis(const dictionary& d);
    void applyMitralSliceAxis();

    // Utilities used by the reader
    void autoFitPlaneFromSurface(vector& nHat, point& centroid) const;
    void maybeFlipNormalUsingValve(vector& nHat, const point& centroid) const;
    void fitTopBottomCapsFromSurface(point& Ct, point& Cb, vector& axisUnit) const;

    // Time/space gains (inline)
    inline scalar timeGain_(scalar t) const
    {
        // Normalized phase x∈[0,1)
        scalar x = 0.0;
        if (mPeriod_ > VSMALL)
        {
            scalar tCycle = std::fmod(t, mPeriod_);
            if (tCycle < 0) tCycle += mPeriod_;
            x = tCycle / mPeriod_;
        }

        // C1 smoothstep
        auto sstep = [](scalar u){ u = max(0.0, min(1.0, u)); return (3*u*u - 2*u*u*u); };

        // helpers with wrap-around support
        const scalar eps = 1e-12;
        auto edgeUp = [&](scalar xx, scalar a, scalar b)->scalar
        {
            if (a <= b)
            {
                if (xx <= a) return 0.0;
                if (xx >= b) return 1.0;
                return sstep((xx - a)/max(eps, b - a));
            }
            if (xx >= a) return sstep((xx - a)/max(eps, 1.0 - a));
            if (xx <= b) return sstep((xx + (1.0 - a))/max(eps, 1.0 - a));
            return 0.0;
        };
        auto edgeDown = [&](scalar xx, scalar a, scalar b)->scalar
        {
            if (a <= b)
            {
                if (xx <= a) return 0.0;
                if (xx >= b) return 1.0;
                return sstep((xx - a)/max(eps, b - a));
            }
            if (xx >= a) return sstep((xx - a)/max(eps, 1.0 - a));
            if (xx <= b) return sstep((xx + (1.0 - a))/max(eps, 1.0 - a));
            return 1.0;
        };

        if (mTimeLaw_ == "twoWindow")
        {
            const scalar rise = edgeUp  (x, mClose0_, mClose1_);
            const scalar fall = edgeDown(x, mOpen0_,  mOpen1_);
            return max(0.0, min(1.0, rise * (1.0 - fall)));
        }
        else if (mTimeLaw_ == "threeWindow")
        {
            const scalar r1 = edgeUp  (x, mClose0_,   mClose1_);
            const scalar fo = edgeDown(x, mOpen0_,    mOpen1_);
            const scalar r2 = edgeUp  (x, mClose2_0_, mClose2_1_);
            return max(0.0, min(1.0, r1*(1.0 - fo) + r2));
        }

        if (mTimeLaw_ == "cos2")
            return 0.5*(1.0 - Foam::cos(2.0*Foam::constant::mathematical::pi*x));
        if (mTimeLaw_ == "smoothstep")
        { scalar s = max(0.0, min(1.0, x)); return (3*s*s - 2*s*s*s); }

        const scalar s = Foam::sin(Foam::constant::mathematical::pi * x);
        return s*s;
    }

    inline scalar spaceGain_(scalar xi) const
    {
        if (xi < 0) return 0.0;
        if (xi <= mNearRamp_)
        {
            const scalar s = clamp(xi / max(VSMALL, mNearRamp_), 0.0, 1.0);
            if (mSpaceLaw_ == "linear")       return s;
            if (mSpaceLaw_ == "smoothstep")   return (3*s*s - 2*s*s*s);   // C1
            /* smootherstep (default) */       return s*s*s*(10.0 - 15.0*s + 6.0*s*s); // C2
        }
        if (xi <= mActiveLen_) return 1.0;
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
