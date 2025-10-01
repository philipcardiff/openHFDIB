#include "createRDeltaT.H"
#include "fvcGrad.H" //<- optional; safe to handle

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field lambda\n" << endl;
volScalarField lambda
(
    IOobject
    (
        "lambda",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh //    dimensionedScalar("zero",dimless,0) // initialising to 0
);


// --- Derived field used by dynamicRefineFvMesh (non-zero near interface)
Info<< "Creating field lambdaRefine (|grad(lambda)|)\n" << endl;
volScalarField lambdaRefine
(
    IOobject
    (
        "lambdaRefine",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mag(fvc::grad(lambda))            // magnitude of gradient -> large only at interface
);
lambdaRefine.correctBoundaryConditions();

// --- cell length scale h ~ (Vcell)^(1/3)
// Build as a volScalarField tied to the mesh, then fill its internal values.
volScalarField hCell
(
    IOobject
    (
        "hCell",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimLength, 0.0)
);

// Fill hCell from cell volumes (works on OF3.x and OF7)
hCell.primitiveFieldRef() = Foam::pow(mesh.V(), 1.0/3.0);
hCell.correctBoundaryConditions();

// --- Normalized, saturated driver in [0,1]: lambdaRefine01 = min(1, h*|grad(lambda)|)
Info<< "Creating field lambdaRefine01 = clamp(h*|grad(lambda)|, 0..1)\n";
volScalarField lambdaRefine01
(
    IOobject("lambdaRefine01",
	     runTime.timeName(),
	     mesh,
             IOobject::NO_READ,
	     IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimless, 0)
);

// clamp to [0,1]
lambdaRefine01 = min(scalar(1), max(scalar(0), hCell*mag(fvc::grad(lambda))));
lambdaRefine01.correctBoundaryConditions();



Info<< "Creating/reading field body\n" << endl;
volScalarField body
(
    IOobject
    (
        "body",                   // <-- exact name used by dynamicMeshDict
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT, // read 0/body if present, else initialise
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0) // default value if not read
);


Info<< "Creating field Ui\n" << endl;
volVectorField Ui
(
    IOobject
    (
        "Ui",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    U
);


volScalarField surface("surf",lambda);

Info<< "Reading field f\n" << endl;
volVectorField f
(
    IOobject
    (
        "f",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
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
