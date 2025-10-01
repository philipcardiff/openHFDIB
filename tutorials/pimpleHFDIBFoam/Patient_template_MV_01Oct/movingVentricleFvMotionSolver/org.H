/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::movingVentricleFvMotionSolver

Group
    grpMeshMotionSolvers

Description
    Mesh motion solver for an fvMesh for a moving ventricle case, where the mesh
    is moved such that the ventricle volume follows a prescribed time curve.

SourceFiles
    movingVentricleFvMotionSolver.C

Author
    Philip Cardiff, UCD.
    Chanikya Valeti, UCD.

\*---------------------------------------------------------------------------*/

#ifndef movingVentricleFvMotionSolver_H
#define movingVentricleFvMotionSolver_H

#include "displacementMotionSolver.H"
#include "fvMotionSolver.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class motionInterpolation;
class motionDiffusivity;

/*---------------------------------------------------------------------------*\
             Class movingVentricleFvMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class movingVentricleFvMotionSolver
:
    public displacementMotionSolver,
    public fvMotionSolver
{
    // Private data

        //- Cell-centre motion field
        mutable volVectorField cellDisplacement_;

        //- Interpolation used to transfer cell displacement to the points
        autoPtr<motionInterpolation> interpolationPtr_;

        //- Diffusivity used to control the motion
        autoPtr<motionDiffusivity> diffusivityPtr_;

        //- Frozen points (that are not on patches). -1 or points that are
        //  fixed to be at points0_ location
        label frozenPointsZone_;

        //- Patch index ventricle moving patch
        const label patchID_;

        //- Minimum volume (assumed to be the starting volume)
        const scalar volumeMin_;

        //- Minimum Y coordinate
        scalar minYPoint_;

        //- Relative tolerance on the volume calculation loop given as a
        //  fraction of the target volume
        const scalar relTol_;

        // Added for Volume integration

        scalarField times_;
        scalarField flowRates_;
        scalarField volumes_;
        bool useFlowIntegration_;

        // Linear interpolation helper
        scalar linearInterp_
        (
            const scalarField& x,
            const scalarField& y,
            const scalar xq
        ) const;


        //- Store the initial patch points
        vectorField initialPatchPoints_;

        //- Store the initial patch point unit normals
        vectorField initialPatchPointNormals_;


    // Private Member Functions

        //- Initialisation steps given a dictionary
        void initialise(const dictionary& dict);

        //- Calculate the domain volume given the patch displacement
        scalar calculateVolume
        (
            const fvMesh& mesh,
            const pointField& initialPoints,
            const labelList& ventricleMeshPoints,
            const vectorField& displacement
        ) const;

        //- No copy construct
        movingVentricleFvMotionSolver
        (
            const movingVentricleFvMotionSolver&
        ) = delete;

        //- No copy assignment
        void operator=(const movingVentricleFvMotionSolver&) = delete;


public:

    //- Runtime type information
    TypeName("movingVentricle");


    // Constructors

        //- Construct from polyMesh and IOdictionary
        movingVentricleFvMotionSolver
        (
            const polyMesh&,
            const IOdictionary&
        );

        //- Construct from components
        movingVentricleFvMotionSolver
        (
            const polyMesh& mesh,
            const IOdictionary& dict,
            const pointVectorField& pointDisplacement,
            const pointIOField& points0
        );


    //- Destructor
    ~movingVentricleFvMotionSolver();


    // Member Functions

        //- Return reference to the diffusivity field
        motionDiffusivity& diffusivity();

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Update topology
        virtual void updateMesh(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
