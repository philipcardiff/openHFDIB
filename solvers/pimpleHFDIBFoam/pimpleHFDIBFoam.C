/* pimpleHFDIBFoam.C — legacy/SDF AMR with unified amrRefine01_outer */

#include "fvCFD.H"
#include "openHFDIB.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "fvcGrad.H"
#include "triSurface.H"
#include "triSurfaceSearch.H"
#include "pointIndexHit.H"
#include "boundBox.H"
#include "IOdictionary.H"
#include "mathematicalConstants.H"
#include <cmath>

// ---------------- util: signed distance (neg solid, pos fluid) -------------
static void computeSignedDistance
(
    const Foam::fvMesh& mesh,
    const Foam::triSurface& surf,
    const Foam::triSurfaceSearch& tss,
    const Foam::volScalarField& body,
    const Foam::scalar outerCut,
    Foam::volScalarField& phiSD
)
{
    using namespace Foam;
    const vector span = 1.1*(mesh.bounds().span() + vector::one*SMALL);

    scalarField& sd = phiSD.primitiveFieldRef();
    const vectorField& C = mesh.C();

    forAll(sd, i)
    {
        const point& p = reinterpret_cast<const point&>(C[i]);
        const pointIndexHit hit = tss.nearest(p, span);
        if (!hit.hit())
        {
            const scalar h = std::cbrt(mesh.V()[i] + SMALL);
            sd[i] = 10.0*h;
            continue;
        }
        const point q = hit.point();
        const scalar d   = mag(q - p);
        const scalar sgn = (body[i] > outerCut ? -1.0 : 1.0);
        sd[i] = sgn*d;
    }
    phiSD.correctBoundaryConditions();
}

// ---------------- util: H_eps, delta_eps, curvature ------------------------
static void buildInterfaceIndicators
(
    const volScalarField& phiSD,
    const volScalarField& hCell,
    const scalar epsFactor,
    volScalarField& Heps,
    volScalarField& Deps,
    volScalarField& kappa
)
{
    using namespace Foam;
    using Foam::constant::mathematical::pi;

    const scalarField& phi = phiSD.primitiveField();
    const scalarField& h   = hCell.primitiveField();
    scalarField& H = Heps.primitiveFieldRef();
    scalarField& D = Deps.primitiveFieldRef();

    forAll(phi, i)
    {
        const scalar e = Foam::max(epsFactor*h[i], 1e-12);
        const scalar x = phi[i];

        if (x <= -e) { H[i] = 0.0; D[i] = 0.0; }
        else if (x >=  e) { H[i] = 1.0; D[i] = 0.0; }
        else
        {
            const scalar r = x/e;
            H[i] = 0.5*(1.0 + r + (1.0/pi)*Foam::sin(pi*r));
            D[i] = 0.5/e * (1.0 + Foam::cos(pi*r));  // ~1/L
        }
    }
    Heps.correctBoundaryConditions();
    Deps.correctBoundaryConditions();

    tmp<volVectorField> gphi = fvc::grad(phiSD);
    const volScalarField magG
    (
        IOobject("magG", phiSD.time().timeName(), phiSD.mesh(),
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mag(gphi())
    );

    volVectorField nhat
    (
        IOobject("nhat", phiSD.time().timeName(), phiSD.mesh(),
                 IOobject::NO_READ, IOobject::NO_WRITE),
        gphi() / max(magG, dimensionedScalar("tiny", dimless, SMALL))
    );

    kappa = fvc::div(nhat);          // ~1/L
    kappa.correctBoundaryConditions();
}

// ==========================================================================

int main(int argc, char *argv[])
{
    #include "postProcess.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createUfIfPresent.H"

    // HFDIB init & first update
    openHFDIB HFDIB(mesh);
    HFDIB.initialize();
    HFDIB.update(body, f);
    body.correctBoundaryConditions();
    // after HFDIB.update(body,f); body.correctBoundaryConditions();
    static bool first = true;
    static volScalarField bodyPrev(body);
    if (!first)
      {
	scalar dSum = gSum(mag(body.primitiveField() - bodyPrev.primitiveField()));
	Info<< "||body - bodyPrev||_1 = " << dSum << nl;
      }
    bodyPrev = body;
    first = false;

    // --- Runtime knobs (controlDict > pimple sub-dict) ---------------------
    const word   amrIndicator      = pimple.dict().lookupOrDefault<word>("amrIndicator", "legacy");
    const scalar amr_epsilonFactor = pimple.dict().lookupOrDefault<scalar>("amr_epsilonFactor", 1.5);
    const scalar amr_outerMaskCut  = pimple.dict().lookupOrDefault<scalar>("amr_outerMaskCut", 0.55);
    const scalar amr_wDelta        = pimple.dict().lookupOrDefault<scalar>("amr_wDelta", 1.0);
    const scalar amr_wCurv         = pimple.dict().lookupOrDefault<scalar>("amr_wCurv", 0.3);
    const scalar amr_kappaClamp    = pimple.dict().lookupOrDefault<scalar>("amr_kappaClamp", 3.0);

    // Read STL for SDF path (once)
    IOdictionary HFDIBDict
    (
        IOobject("HFDIBDict", runTime.constant(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE)
    );
    wordList bodyNames(HFDIBDict.lookup("bodyList"));
    const word bodyName = bodyNames[0];
    const dictionary& bodyDict = HFDIBDict.subDict(bodyName);
    Foam::fileName rawName = Foam::fileName(bodyDict.lookup("fileName"));
    Foam::fileName stlPath = (rawName.isAbsolute() ? rawName : Foam::fileName("constant/triSurface")/rawName);
    Foam::triSurface surf(stlPath);
    Foam::triSurfaceSearch tss(surf);

    // Helper: build AMR indicators (legacy or SDF) and update amrRefine01_outer
    auto buildAMR = [&]()
    {
        // Refresh fluid-only mask with chosen cut
        {
            scalarField& m = outerMask.primitiveFieldRef();
            const scalarField& b = body.primitiveFieldRef();
            forAll(m, i) m[i] = (b[i] < amr_outerMaskCut ? 1.0 : 0.0);
            outerMask.correctBoundaryConditions();
        }

        if (amrIndicator == "sdf")
        {
            // --- SDF-based (robust)
            computeSignedDistance(mesh, surf, tss, body, amr_outerMaskCut, phiSD);
            buildInterfaceIndicators(phiSD, hCell, amr_epsilonFactor, Heps, Deps, kappa);

            etaRefine01 =
                min( scalar(1),
                     max( scalar(0),
                          amr_wDelta * (Deps * (amr_epsilonFactor*hCell))    // δ_ε * ε → [-]
                        + amr_wCurv  * min( amr_kappaClamp, mag(kappa)*hCell ) // κ*h → [-]
                     ) );
            etaRefine01.correctBoundaryConditions();

            etaRefine01_outer = etaRefine01 * outerMask;
            etaRefine01_outer.correctBoundaryConditions();

            amrRefine01_outer = etaRefine01_outer;   // <-- SWITCH target (SDF)
            amrRefine01_outer.correctBoundaryConditions();
        }
        else
        {
            // --- Legacy: lambda <- body; smooth; |grad|; scale; clamp; mask
            lambda = body;
            lambda.correctBoundaryConditions();

            lambdaSmooth = lambda;
            for (label i=0; i<2; ++i) { fvc::smooth(lambdaSmooth, 0.8); }
            lambdaSmooth = min(scalar(1), max(scalar(0), lambdaSmooth));
            lambdaSmooth.correctBoundaryConditions();

            lambdaRefine  = mag(fvc::grad(lambdaSmooth));   // 1/L
            lambdaRefine.correctBoundaryConditions();

            lambdaRefine01 = min(scalar(1), max(scalar(0), hCell*lambdaRefine)); // [-]
            lambdaRefine01.correctBoundaryConditions();

            lambdaRefine01_outer = lambdaRefine01 * outerMask; // [-]
            lambdaRefine01_outer.correctBoundaryConditions();

            amrRefine01_outer = lambdaRefine01_outer;       // <-- SWITCH target (legacy)
            amrRefine01_outer.correctBoundaryConditions();
        }
    };

    // Initial AMR build on the start mesh
    buildAMR();
    turbulence->validate();

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "setDeltaT.H"
        }

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Keep IB fields current on the present mesh
        HFDIB.update(body, f);
        body.correctBoundaryConditions();

        // Rebuild AMR driver for this step
        buildAMR();

        // --- PIMPLE loop with dynamic mesh
        while (pimple.loop())
        {
            HFDIB.update(body, f);
            body.correctBoundaryConditions();

            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                // Dynamic mesh update (driven by amrRefine01_outer in dynamicMeshDict)
                mesh.update();

                if (mesh.changing())
                {
                    MRF.update();

                    if (correctPhi)
                    {
                        // absolute flux from mapped surface velocity
                        phi = mesh.Sf() & Uf();
                        #include "CorrectPhi.H"
                        // make flux relative to mesh motion
                        fvc::makeRelative(phi, U);
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }

                    // Rebuild IB & AMR on the NEW mesh (volumes changed)
                    HFDIB.update(body, f);
                    body.correctBoundaryConditions();

                    // recompute hCell on the new mesh
                    {
                        scalarField& h = hCell.primitiveFieldRef();
                        const scalarField& V = mesh.V().field();
                        forAll(h, i) { h[i] = std::cbrt( max(V[i], VSMALL) ); }
                    }
                    hCell.correctBoundaryConditions();

                    buildAMR();
                }
            }

            #include "UEqn.H"

            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl << endl;
    }

    Info<< "End\n" << endl;
    return 0;
}
