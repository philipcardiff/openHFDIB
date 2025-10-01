//* createFields.H — unified AMR driver (legacy or SDF) via amrRefine01_outer *//

#include "createRDeltaT.H"
#include "fvcGrad.H"          // for fvc::grad
#include "fvcSmooth.H"        // for fvc::smooth
#include "IOdictionary.H"
#include <cmath>

//---------------------------------------------------------------------------//
// Primary fields
//---------------------------------------------------------------------------//

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject("p", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject("U", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh
);

//---------------------------------------------------------------------------//
// Base IB indicator from STL (0 fluid, 1 solid, ~0.5 interface)
//---------------------------------------------------------------------------//

Info<< "Creating/reading field body\n" << endl;
volScalarField body
(
    IOobject("body", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

// Fluid-only mask (1 in fluid, 0 in/near solid). Slight bias toward fluid.
volScalarField outerMask
(
    IOobject("outerMask", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

// Build outerMask = 1 if body < 0.55 else 0 (solver may overwrite each step)
{
    scalarField& m = outerMask.primitiveFieldRef();
    const scalarField& b = body.primitiveFieldRef();
    forAll(m, i) { m[i] = (b[i] < 0.55 ? 1.0 : 0.0); }
}
outerMask.correctBoundaryConditions();

//---------------------------------------------------------------------------//
// Legacy indicator scaffolding: lambda = body → smooth → clamp
//---------------------------------------------------------------------------//

Info<< "Reading/creating field lambda (legacy indicator)\n" << endl;
volScalarField lambda
(
    IOobject("lambda", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

volScalarField lambdaSmooth
(
    IOobject("lambdaSmooth", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

// initialise lambda from body and smooth
lambda = body;
lambda.correctBoundaryConditions();

lambdaSmooth = lambda;
for (label i=0; i<2; ++i) { fvc::smooth(lambdaSmooth, 0.8); }  // 2 passes is robust
lambdaSmooth = min(scalar(1), max(scalar(0), lambdaSmooth));
lambdaSmooth.correctBoundaryConditions();

//---------------------------------------------------------------------------//
// Legacy gradients & scaling
//---------------------------------------------------------------------------//

Info<< "Creating field lambdaRefine = |grad(lambdaSmooth)|\n" << endl;
volScalarField lambdaRefine
(
    IOobject("lambdaRefine", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mag(fvc::grad(lambdaSmooth))                    // dimensions: 1/length
);
// NOTE: expression constructor ⇒ no extra 'mesh' argument here
lambdaRefine.correctBoundaryConditions();

// Cell size h ≈ cbrt(volume). Recomputed again after AMR in the solver.
volScalarField hCell
(
    IOobject("hCell", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("zero", dimLength, 0.0)
);
{
    scalarField& h = hCell.primitiveFieldRef();
    const scalarField& V = mesh.V().field();   // raw numeric volumes
    forAll(h, i) { h[i] = std::cbrt(Foam::max(V[i], VSMALL)); }
}
hCell.correctBoundaryConditions();

// Dimensionless 0..1 legacy indicator: clamp( h * |∇lambdaSmooth| )
Info<< "Creating field lambdaRefine01 = clamp(h*|grad(lambdaSmooth)|,0..1)\n" << endl;
volScalarField lambdaRefine01
(
    IOobject("lambdaRefine01", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    min(scalar(1), max(scalar(0), hCell*lambdaRefine))       // dimensions: L*(1/L) → 1
);
// NOTE: expression constructor ⇒ no extra 'mesh' argument here
lambdaRefine01.correctBoundaryConditions();

// Fluid-only legacy driver
volScalarField lambdaRefine01_outer
(
    IOobject("lambdaRefine01_outer", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    lambdaRefine01 * outerMask
);
lambdaRefine01_outer.correctBoundaryConditions();

//---------------------------------------------------------------------------//
// SDF-based AMR scaffolding
//---------------------------------------------------------------------------//

Info<< "Creating SDF-based AMR fields (phiSD, Heps, Deps, kappa, etaRefine01, etaRefine01_outer)\n" << endl;

volScalarField phiSD
(
    IOobject("phiSD", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimLength, 0.0)                // signed distance [L]
);

volScalarField Heps
(
    IOobject("Heps", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)                  // Heaviside H_ε [-]
);

volScalarField Deps
(
    IOobject("Deps", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimless/dimLength, 0.0)        // δ_ε ~ 1/L
);

volScalarField kappa
(
    IOobject("kappa", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimless/dimLength, 0.0)        // curvature ~ 1/L
);

volScalarField etaRefine01
(
    IOobject("etaRefine01", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)                  // [-]
);

volScalarField etaRefine01_outer
(
    IOobject("etaRefine01_outer", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)                  // [-]
);

//---------------------------------------------------------------------------//
// Unified AMR driver (switch target; solver overwrites each step)
//---------------------------------------------------------------------------//

Info<< "Creating field amrRefine01_outer (switchable target for AMR)\n" << endl;
volScalarField amrRefine01_outer
(
    IOobject("amrRefine01_outer", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    lambdaRefine01_outer                                      // default to legacy
);
amrRefine01_outer.correctBoundaryConditions();

//---------------------------------------------------------------------------//
// Misc standard fields
//---------------------------------------------------------------------------//

Info<< "Creating field Ui\n" << endl;
volVectorField Ui
(
    IOobject("Ui", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
    U
);

Info<< "Reading field f\n" << endl;
volVectorField f
(
    IOobject("f", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh
);

#include "createPhi.H"

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

#include "createMRF.H"
#include "createFvOptions.H"
