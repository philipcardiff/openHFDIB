/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "movingVentricleFvMotionSolver.H"
#include "motionInterpolation.H"
#include "motionDiffusivity.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "surfaceInterpolate.H"
#include "fvcLaplacian.H"
#include "mapPolyMesh.H"
#include "fvOptions.H"
#include "IFstream.H"
#include "PrimitivePatch.H"
#include "PatchTools.H"
#include <cmath>
using std::isfinite;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(movingVentricleFvMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        movingVentricleFvMotionSolver,
        dictionary
    );

    addToRunTimeSelectionTable
    (
        displacementMotionSolver,
        movingVentricleFvMotionSolver,
        displacement
    );
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::movingVentricleFvMotionSolver::initialise(const dictionary& dict)
{
    if (patchID_ == -1)
    {
        FatalErrorInFunction
            << "Did not find the patch = " << word(dict.lookup("patchName"))
            << endl;
    }

    Info<< typeName<< ": minimum volume = " << volumeMin_ << endl;


    // Read time vs flow rate series
    interpolationTable<scalar> flowTable
    (
        motionSolver::coeffDict().subDict("timeVsFlowSeries")
    );

    // Extract times and flow values
    times_.setSize(flowTable.size());
    flowRates_.setSize(flowTable.size());

    {
        label i = 0;
        forAllConstIter(interpolationTable<scalar>, flowTable, iter)
        {
            const Tuple2<scalar, scalar>& tp = *iter;
            times_[i] = tp.first();
            flowRates_[i] = tp.second();
            ++i;
        }
    }

    // Integrate the flow rate  using trapezoidal rule to get the volume versus
    // time
    volumes_.setSize(times_.size(), 0.0);
    for (label i = 1; i < times_.size(); ++i)
    {
        const scalar dt = times_[i] - times_[i-1];
        const scalar avgFlow = 0.5*(flowRates_[i] + flowRates_[i - 1]);
        volumes_[i] = volumes_[i - 1] + avgFlow*dt;
    }

    // Offset by minimum volume
    forAll(volumes_, i)
    {
        volumes_[i] += volumeMin_;
    }

    Info<< typeName << ": using timeVsFlowSeries with "
        << times_.size() << " samples, volumeMin = "
        << volumeMin_ << endl;

    // Set min Y point
    minYPoint_ =
        gMin
        (
            fvMesh_.boundaryMesh()[patchID_].localPoints().component(vector::Y)
        );

    // Set the initial patch points
    initialPatchPoints_ = fvMesh_.boundaryMesh()[patchID_].localPoints();

    // Set the initial point normals
    // Ensure we use parallel-consistent normals
    initialPatchPointNormals_ =
        PatchTools().pointNormals(fvMesh_, fvMesh_.boundaryMesh()[patchID_]);
}


Foam::scalar Foam::movingVentricleFvMotionSolver::linearInterp
(
    const scalarField& x, const scalarField& y, const scalar xq
) const
{
    if (xq <= x.first())
    {
        return y.first();
    }

    if (xq >= x.last())
    {
        return y.last();
    }

    label i = 1;
    for (; i < x.size(); ++i)
    {
        if (x[i] >= xq) break;
    }

    const scalar x0 = x[i -1], x1 = x[i];
    const scalar y0 = y[i -1], y1 = y[i];
    const scalar w  = (xq - x0)/max(VSMALL, x1 - x0);

    return y0 + w*(y1 - y0);
}


Foam::scalar Foam::movingVentricleFvMotionSolver::calculateVolume
(
    const fvMesh& mesh,
    const pointField& initialPoints,
    const labelList& ventricleMeshPoints,
    const vectorField& displacement
) const
{
    // We are going to do something a bit hacky: we will use const_cast to move
    // the mesh by the the given displacement field, then calculate the new
    // volume, and finally reset the mesh before returning the volume

    // Take a copy of the old time mesh points
    const vectorField oldPoints(mesh.points());

    // Calculate the new points
    pointField newPoints(initialPoints);
    forAll(ventricleMeshPoints, pI)
    {
        const label pointID = ventricleMeshPoints[pI];
        newPoints[pointID] += displacement[pI];
    }

    // Move the mesh
    const_cast<fvMesh&>(mesh).movePoints(newPoints);

    // Calculate the volume
    const scalar totalVolume = gSum(mesh.cellVolumes());

    // Reset the mesh
    const_cast<fvMesh&>(mesh).movePoints(oldPoints);

    return totalVolume;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::movingVentricleFvMotionSolver::movingVentricleFvMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    displacementMotionSolver(mesh, dict, typeName),
    fvMotionSolver(mesh),
    cellDisplacement_
    (
        IOobject
        (
            "cellDisplacement",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvMesh_,
        dimensionedVector(pointDisplacement_.dimensions(), Zero),
        cellMotionBoundaryTypes<vector>(pointDisplacement_.boundaryField())
    ),
    interpolationPtr_
    (
        coeffDict().found("interpolation")
      ? motionInterpolation::New(fvMesh_, coeffDict().lookup("interpolation"))
      : motionInterpolation::New(fvMesh_)
    ),
    diffusivityPtr_
    (
        motionDiffusivity::New(fvMesh_, coeffDict().lookup("diffusivity"))
    ),
    frozenPointsZone_
    (
        coeffDict().found("frozenPointsZone")
      ? fvMesh_.pointZones().findZoneID
        (
            coeffDict().get<word>("frozenPointsZone")
        )
      : -1
    ),
    patchID_
    (
        mesh.boundaryMesh().findPatchID
        (
            word(motionSolver::coeffDict().lookup("patchName"))
        )
    ),
    volumeMin_(gSum(mesh.cellVolumes())),
    minYPoint_(0.0),
    useNewtonRaphson_
    (
        motionSolver::coeffDict().lookupOrDefault<Switch>
        (
            "useNewtonRaphson", true
        )
    ),
    relTol_
    (
        motionSolver::coeffDict().lookupOrDefault<scalar>("relTol", 1e-6)
    ),
    initialPatchPoints_(),
    initialPatchPointNormals_()
{
    initialise(dict);
}


Foam::movingVentricleFvMotionSolver::movingVentricleFvMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict,
    const pointVectorField& pointDisplacement,
    const pointIOField& points0
)
:
    displacementMotionSolver(mesh, dict, pointDisplacement, points0, typeName),
    fvMotionSolver(mesh),
    cellDisplacement_
    (
        IOobject
        (
            "cellDisplacement",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvMesh_,
        dimensionedVector(pointDisplacement_.dimensions(), Zero),
        cellMotionBoundaryTypes<vector>(pointDisplacement_.boundaryField())
    ),
    interpolationPtr_
    (
        coeffDict().found("interpolation")
      ? motionInterpolation::New(fvMesh_, coeffDict().lookup("interpolation"))
      : motionInterpolation::New(fvMesh_)
    ),
    diffusivityPtr_
    (
        motionDiffusivity::New(fvMesh_, coeffDict().lookup("diffusivity"))
    ),
    frozenPointsZone_
    (
        coeffDict().found("frozenPointsZone")
      ? fvMesh_.pointZones().findZoneID
        (
            coeffDict().get<word>("frozenPointsZone")
        )
      : -1
    ),
    patchID_
    (
        mesh.boundaryMesh().findPatchID
        (
            word(motionSolver::coeffDict().lookup("patchName"))
        )
    ),
    volumeMin_(gSum(mesh.cellVolumes())),
    minYPoint_(0.0),
    useNewtonRaphson_
    (
        motionSolver::coeffDict().lookupOrDefault<Switch>
        (
            "useNewtonRaphson", true
        )
    ),
    relTol_
    (
        motionSolver::coeffDict().lookupOrDefault<scalar>("relTol", 1e-3)
    ),
    initialPatchPoints_(),
    initialPatchPointNormals_()
{
    initialise(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::movingVentricleFvMotionSolver::~movingVentricleFvMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::motionDiffusivity&
Foam::movingVentricleFvMotionSolver::diffusivity()
{
    if (!diffusivityPtr_)
    {
        diffusivityPtr_ = motionDiffusivity::New
        (
            fvMesh_,
            coeffDict().lookup("diffusivity")
        );
    }

    return *diffusivityPtr_;
}


Foam::tmp<Foam::pointField>
Foam::movingVentricleFvMotionSolver::curPoints() const
{
    interpolationPtr_->interpolate
    (
        cellDisplacement_,
        pointDisplacement_
    );

    tmp<pointField> tcurPoints
    (
        points0() + pointDisplacement_.primitiveField()
    );
    pointField& curPoints = tcurPoints.ref();

    // Implement frozen points
    if (frozenPointsZone_ != -1)
    {
        const pointZone& pz = fvMesh_.pointZones()[frozenPointsZone_];

        forAll(pz, i)
        {
            curPoints[pz[i]] = points0()[pz[i]];
        }
    }

    twoDCorrectPoints(curPoints);

    return tcurPoints;
}


void Foam::movingVentricleFvMotionSolver::solve()
{
    // The points have moved so before interpolation update
    // the motionSolver accordingly
    movePoints(fvMesh_.points());

    const Time& runTime = fvMesh_.time();

    // Take a reference to the patch and its point indices
    const polyPatch& ventriclePatch = fvMesh_.boundaryMesh()[patchID_];
    const labelList& ventricleMeshPoints = ventriclePatch.meshPoints();

    // Lookup the current volume
    const scalar targetVolume = linearInterp(times_, volumes_, runTime.value());

    // Determine the displacement scaling factor to achieve the target volume
    if (useNewtonRaphson_)
    {
        // Initialise the patch point displacement field
        vectorField patchPointDisp(ventriclePatch.nPoints(), vector::zero);

        // Geometry scale (so the target is reachable)
        const scalar kScale = 1.0;

        // Apply scale 'a', enforce BCs, return volume
        auto evalVolume = [&](scalar a) -> scalar
        {
            vectorField patchPointDisp
            (
                initialPatchPointNormals_
               *Foam::max(0.0, initialPatchPoints_.component(vector::Y)/minYPoint_)
               *kScale*a
            );

            // Push to boundary, enforce continuity across procs
            pointDisplacement_.boundaryFieldRef()[patchID_] == patchPointDisp;
            pointDisplacement_.correctBoundaryConditions();
            patchPointDisp =
                pointDisplacement_.boundaryFieldRef()
                [
                    patchID_
                ].patchInternalField();

            return
                calculateVolume
                (
                    fvMesh_, points0(), ventricleMeshPoints, patchPointDisp
                );
        };

        // Newton-Raphosn using forward finite differencing to approximate the
        // deriviative
        scalar a = 0.0;
        scalar V = evalVolume(a);

        const scalar tol = relTol_*targetVolume;
        const label  maxIter = 30;
        label iter = 0;
        scalar res =  mag(V - targetVolume);

        Info<< "Find the displacement scale factor using the Newton-Raphson"
            << "method" << nl
            << "    Target volume = " << targetVolume << endl;

        while (res > tol && ++iter < maxIter)
        {
            // Forward-difference slope: dV/da ≈ (V(a+h) - V(a)) / h
            const scalar h = Foam::max(1e-8*mag(a), 1e-8*kScale);
            const scalar Vp = evalVolume(a + h);
            const scalar dV = (Vp - V) / h;

            if (mag(dV) < VSMALL)
            {
                Info<< "  Newton: tiny slope, stopping.\n";
                break;
            }

            scalar aNew = a - (V - targetVolume)/dV;

            // Light damping: ensure residual decreases; at most 5 halvings
            scalar Vnew = evalVolume(aNew);
            label ls = 0;
            while
            (
                mag(Vnew - targetVolume) > 0.9*mag(V - targetVolume)
             && ls++ < 5
            )
            {
                aNew = 0.5*(a + aNew);
                Vnew = evalVolume(aNew);
            }

            a = aNew;
            V = Vnew;

            // Update the residual
            res =  mag(V - targetVolume);

            // Print info
            Info<< "    " << iter
                << ": scale factor = " << a
                << ", volume = " << V
                << ", |res| = " << res << nl;

            if (iter == maxIter)
            {
                FatalErrorInFunction
                    << "Max iterations reached in the Newton-Raphson loop!"
                    << exit(FatalError);
            }
        }

        // Final: keep the converged displacement in patchPointDisp
        patchPointDisp =
            initialPatchPointNormals_
           *Foam::max(0.0, initialPatchPoints_.component(vector::Y)/minYPoint_)
           *kScale*a;

        // Push to boundary, enforce continuity across procs
        pointDisplacement_.boundaryFieldRef()[patchID_] == patchPointDisp;
        pointDisplacement_.correctBoundaryConditions();
        patchPointDisp =
            pointDisplacement_.boundaryFieldRef()[patchID_].patchInternalField();

        const scalar totalVolume =
            calculateVolume
            (
                fvMesh_, points0(), ventricleMeshPoints, patchPointDisp
            );

        Info<< "    Converged: scale factor = " << a
            << ", volume = " << totalVolume
            << ", |res| = " << res << endl;
    }
    else // Bracketed bisection meth0d
    {
        // Lookup the old target volume
        const scalar prevVolume =
        linearInterp(times_, volumes_, runTime.value() - runTime.deltaTValue());

        // Choose the monotonic direction for initial bracket
        scalar aLow = 0;
        scalar aHigh = 0;
        if (prevVolume < targetVolume)
        {
            aLow = 0.0;
            aHigh =  1.0;
        }
        else if (prevVolume > targetVolume)
        {
            aLow = -1.0;
            aHigh = 0.0;
        }
        else
        {
            aLow = 0.0;
            aHigh = 0.0;
        }

        // Helper to update the patchPointDisp field given a factor 'a'
        auto updatePatchPointDisp = [&](vectorField& patchPointDisp, scalar a) -> void
        {
            // aHighMax controls the max volume guess
            // If it is too small, then the solution will never be reached
            // If too big, then many iterations may be required
            const scalar aHighMax = 1.0;

            // Build patchPointDisp for this 'a'
            patchPointDisp =
                initialPatchPointNormals_
               *Foam::max(0.0, initialPatchPoints_.component(vector::Y)/minYPoint_)
               *aHighMax*a;
        };

        // Helper to apply factor 'a' and compute volume
        auto evalVolume = [&](scalar a) -> scalar
        {
            // Build patchPointDisp for this 'a'
            vectorField patchPointDisp
            (
                initialPatchPointNormals_.size(), vector::zero
            );
            updatePatchPointDisp(patchPointDisp, a);

            // Push to boundary, enforce continuity across procs
            pointDisplacement_.boundaryFieldRef()[patchID_] == patchPointDisp;
            pointDisplacement_.correctBoundaryConditions();
            patchPointDisp =
                pointDisplacement_.boundaryFieldRef()
                [
                    patchID_
                ].patchInternalField();

            // Compute volume
            return
                calculateVolume
                (
                    fvMesh_, points0(), ventricleMeshPoints, patchPointDisp
                );
        };

        // Evaluate at a = 0 once (cheap baseline and often close to target)
        scalar Vlow = evalVolume(aLow);
        scalar V0 = Vlow;

        Info<< "Calculating the " << ventriclePatch.name() << " patch motion:" << nl
            << "    targetVolume = " << targetVolume << nl
            << "    V(0)         = " << V0 << endl;

        // If already close, stop
        scalar aTry = 0.0;
        if (mag(targetVolume - V0) <= relTol_*targetVolume)
        {
            Info<< "    Already within tolerance at a = 0" << endl;
            // Ensure the zero displacement is set (already done by evalVolume(aLow))
            // Exit early
        }
        else
        {
            // Ensure we have a proper bracket [aLow, aHigh]
            // Expand aHigh (or aLow) geometrically until it brackets the target
            const scalar expandFactor = 2.0;
            const scalar maxAbsA = 10.0;   // hard safety cap on factor magnitude
            const int maxExpand = 20;

            // Pick an initial opposite end if degenerate
            if (aLow == aHigh)
            {
                aHigh = (targetVolume > Vlow ? 1.0 : -1.0);
            }

            scalar Vhigh = evalVolume(aHigh);
            int expandI = 0;
            auto bracketed = [&](scalar Va, scalar Vb)
            {
                return ( (Va - targetVolume) * (Vb - targetVolume) <= 0 );
            };

            while (!bracketed(Vlow, Vhigh) && ++expandI <= maxExpand)
            {
                // Expand toward the side closer to target
                if
                (
                    Foam::mag(Vhigh - targetVolume) < Foam::mag(Vlow - targetVolume)
                )
                {
                    aHigh =
                        Foam::sign(aHigh)
                       *Foam::min
                        (
                            maxAbsA, expandFactor*Foam::max(1e-6, Foam::mag(aHigh))
                        );
                    Vhigh = evalVolume(aHigh);
                }
                else
                {
                    aLow =
                        Foam::sign(aLow)
                       *Foam::min
                        (
                            maxAbsA, expandFactor*Foam::max(1e-6, Foam::mag(aLow))
                        );
                    Vlow = evalVolume(aLow);
                }
            }

            if (!bracketed(Vlow, Vhigh))
            {
                WarningInFunction
                    << "Failed to bracket target volume after expansion. "
                    << "Proceeding with cautious bisection around a=0." << nl;

                // Fall back to symmetric bracket around 0
                aLow  = -1.0;
                Vlow  = evalVolume(aLow);
                aHigh =  1.0;
                Vhigh = evalVolume(aHigh);
            }

            // Safeguarded secant/bisection + optional ratio step
            const label maxIter = 50;
            const scalar tol = relTol_*targetVolume;
            const scalar gamma = 1.0;      // ratio exponent (0.5–1.0 is fine)
            const scalar stepLimit = 3.0;  // limit change per iteration (×)
            label iter = 0;

            while (++iter <= maxIter)
            {
                // Midpoint (always inside)
                const scalar amid = 0.5*(aLow + aHigh);
                scalar Vmid = evalVolume(amid);
                Info<< "    " << iter << ": a∈[" << aLow << ", " << aHigh
                    << "], amid = " << amid << ", V = " << Vmid << endl;

                if (mag(Vmid - targetVolume) <= tol)
                {
                    Info<< "    Converged: V = " << Vmid << " at a = " << amid << endl;
                    break;
                }

                // --- Secant proposal ---
                scalar aSec = amid;
                {
                    const scalar dV = (Vhigh - Vlow);
                    if (Foam::mag(dV) > SMALL)
                    {
                        // Linear interpolation for root between (aLow,Vlow) and
                        // (aHigh,Vhigh)
                        aSec = aHigh - (Vhigh - targetVolume)*(aHigh - aLow)/dV;
                    }
                }

                // --- Ratio proposal (multiplicative) ---
                scalar aRatio = amid;
                {
                    // Use current Vmid; scale a by (target/V)^gamma
                    const scalar s =
                        Foam::max
                        (
                            0.1,
                            Foam::min
                            (
                                10.0,
                                Foam::pow
                                (
                                    targetVolume/Foam::max(SMALL, Vmid),
                                    gamma
                                )
                            )
                        );
                    const scalar deltaLimit =
                        Foam::max(1.0/stepLimit, Foam::min(stepLimit, s));
                    aRatio = amid*deltaLimit;
                }

                // Blend proposals, then **safeguard into the bracket**
                aTry = 0.7*aSec + 0.3*aRatio;

                // If secant went out of bracket or NAN/INF, fall back to bisection
                auto inBracket = [&](scalar a)
                {
                    return (Foam::min(aLow, aHigh) <= a && a <= Foam::max(aLow, aHigh));
                };

                if (!std::isfinite(aTry) || !inBracket(aTry))
                {
                    aTry = amid; // bisection
                }

                // Evaluate and update bracket
                const scalar Vtry = evalVolume(aTry);
                if ((Vlow - targetVolume) * (Vtry - targetVolume) <= 0)
                {
                    aHigh = aTry;
                    Vhigh = Vtry;
                }
                else
                {
                    aLow  = aTry;
                    Vlow  = Vtry;
                }

                if (Foam::mag(aHigh - aLow) <= 1e-6) // tiny interval in 'a'
                {
                    Info<< "    Bracket collapsed in 'a'." << endl;
                    break;
                }
            }
        }

        // Assign the patchPointDisp field to the pointDisplacement patch boundary
        // condition
        vectorField patchPointDisp
        (
            initialPatchPointNormals_.size(), vector::zero
        );
        updatePatchPointDisp(patchPointDisp, aTry);
        pointDisplacement_.boundaryFieldRef()[patchID_] == patchPointDisp;
        pointDisplacement_.correctBoundaryConditions();
    }

    // Update the diffusivity field
    diffusivity().correct();

    // Update the pointDisplacement boundary conditions: these will be used
    // by the cellDisplacement boundary conditions
    pointDisplacement_.boundaryFieldRef().updateCoeffs();

    // Solve for cellDisplacement
    Info<< "Solving the mesh motion for cellDisplacement" << endl;
    fvVectorMatrix DEqn
    (
        fvm::laplacian
        (
            dimensionedScalar("viscosity", dimViscosity, 1.0)
           *diffusivity().operator()(),
            cellDisplacement_,
            "laplacian(diffusivity,cellDisplacement)"
        )
    );

    DEqn.solveSegregatedOrCoupled();
}


void Foam::movingVentricleFvMotionSolver::updateMesh
(
    const mapPolyMesh& mpm
)
{
    displacementMotionSolver::updateMesh(mpm);
}



// ************************************************************************* //
