/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.co
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

// Added here 

useFlowIntegration_ = false;

if (motionSolver::coeffDict().found("timeVsFlowSeries"))
{
    const dictionary& flowDict =
        motionSolver::coeffDict().subDict("timeVsFlowSeries");

    // Load from file
    interpolationTable<scalar> flowTable(flowDict);

    // Extract times and flow values
    times_.setSize(flowTable.size());
    flowRates_.setSize(flowTable.size());

    label i = 0;
    forAllConstIter(interpolationTable<scalar>, flowTable, iter)
    {
        label i = 0;
        for
        (
             interpolationTable<scalar>::const_iterator iter = flowTable.begin();
             iter != flowTable.end();
             ++iter, ++i
        )
    {
    const Tuple2<scalar, scalar>& tp = *iter;
    times_[i]     = tp.first();   // time
    flowRates_[i] = tp.second();  // flow
    }

        ++i;
    }

    // Integrate using trapezoidal rule
    volumes_.setSize(times_.size(), 0.0);
    for (label i = 1; i < times_.size(); ++i)
    {
        scalar dt = times_[i] - times_[i-1];
        scalar avgFlow = 0.5*(flowRates_[i] + flowRates_[i-1]);
        volumes_[i] = volumes_[i-1] + avgFlow*dt;
    }

    // Offset by minimum volume
    forAll(volumes_, i)
    {
        volumes_[i] += volumeMin_;
    }

    useFlowIntegration_ = true;

    Info<< typeName << ": using timeVsFlowSeries with "
        << times_.size() << " samples, volumeMin = "
        << volumeMin_ << endl;
}

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
    relTol_
    (
        motionSolver::coeffDict().lookupOrDefault<scalar>("relTol", 1e-3)
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

    //Interpolating the volume
  
    scalar targetVolume = linearInterp_(times_, volumes_, runTime.value());

    // Initialise the previous volume as the minimum volume
    scalar prevVolume = volumeMin_;

    //  Is the volume increasing, decreasing or equal?
    scalar displacementFactorLow = 0;
    scalar displacementFactorMid = 0;
    scalar displacementFactorHigh = 0;
    if (prevVolume < targetVolume)
    {
        displacementFactorLow = 0;
        displacementFactorHigh = 1;
    }
    else if (prevVolume > targetVolume)
    {
        displacementFactorLow = -1;
        displacementFactorHigh = 0;
    }
    else
    {
        displacementFactorLow = 0;
        displacementFactorHigh = 0;
    }

    // Allocate the displacement field for each patch point
    vectorField patchPointDisp(ventriclePatch.nPoints(), vector::zero);

    // Initialise the total volume
    scalar totalVolume = 0;

    // Iteratively move the boundary patch until the volume matches the target
    // volume
    label iteration = 0;
    const label maxIter = 100;
    Info<< "Calculating the " << ventriclePatch.name() << " patch motion:" << nl
        << "    targetVolume = " << targetVolume << endl;
    while
    (
        mag(targetVolume - totalVolume) > relTol_*targetVolume
     && ++iteration < maxIter
    )
    {
        // Calculate the patchPointDisp for each vertex using midpoint of the
        // interval
        displacementFactorMid =
            (displacementFactorLow + displacementFactorHigh)/2.0;

        // Calculate patchPointDisp field
        patchPointDisp =
            initialPatchPointNormals_
           *Foam::max
            (
                0.0,
                initialPatchPoints_.component(vector::Y)/minYPoint_
            )
           *1e-2*displacementFactorMid;

        // We will copy the patch patchPointDisp field to the pointDisplacement
        // field to enforce continuity between processors
        // Note: the "==" overwrites the patch value field
        pointDisplacement_.boundaryFieldRef()[patchID_] == patchPointDisp;
        pointDisplacement_.correctBoundaryConditions();
        patchPointDisp =
            pointDisplacement_.boundaryFieldRef()
            [
                patchID_
            ].patchInternalField();

        // Calculate the volume of the ventricle given on the latest
        // patchPointDisp field
        totalVolume =
            calculateVolume
            (
                fvMesh_, points0(), ventricleMeshPoints, patchPointDisp
            );

        Info<< "    " << iteration << ": totalVolume = " << totalVolume << endl;

        if (totalVolume > targetVolume)
        {
            displacementFactorHigh = displacementFactorMid;
        }
        else
        {
            displacementFactorLow = displacementFactorMid;
        }
    }

    if (iteration == (maxIter - 1))
    {
        FatalErrorInFunction
            << "Maximum number of iterations reached!" << exit(FatalError);
    }

    // Assign the patchPointDisp field to the pointDisplacement patch boundary
    // condition
    pointDisplacement_.boundaryFieldRef()[patchID_] == patchPointDisp;
    pointDisplacement_.correctBoundaryConditions();

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

// Added for Interpolation purpose
Foam::scalar Foam::movingVentricleFvMotionSolver::linearInterp_
(
    const scalarField& x, const scalarField& y, const scalar xq
) const
{
    if (xq <= x.first()) return y.first();
    if (xq >= x.last())  return y.last();

    label i = 1;
    for (; i < x.size(); ++i)
    {
        if (x[i] >= xq) break;
    }

    scalar x0 = x[i-1], x1 = x[i];
    scalar y0 = y[i-1], y1 = y[i];
    scalar w  = (xq - x0)/max(VSMALL, x1 - x0);

    return y0 + w*(y1 - y0);
}


// ************************************************************************* //
