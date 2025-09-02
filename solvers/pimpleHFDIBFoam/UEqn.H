// Update the Ui and surface fields
{
    vectorField& UiI = Ui.primitiveFieldRef();;
    scalarField& surfaceI = surface.primitiveFieldRef();;
    const scalarField& lambdaI = lambda.primitiveFieldRef();;
    forAll(UiI, cellI)
    {
        UiI[cellI].component(0) = 0.0; // Why do they set Ui.x to zero?!

        if (lambdaI[cellI] > 0.001)
        {
            surfaceI[cellI] = 1.0;
        }
        else
        {
            surfaceI[cellI] = 0.0;
        }
    }

    // Processor boundaries need to be updated
    Ui.correctBoundaryConditions();
    surface.correctBoundaryConditions();
}

// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(U) + fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevReff(U)
 ==
    fvOptions(U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

// Update the immersed force
volScalarField rAU(1.0/UEqn.A());
HFDIB.interpolateIB(U, Ui, lambda);
f = surface*((1.0/rAU)*Ui - UEqn.H() + fvc::grad(p));

//if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p) + f);

    fvOptions.correct(U);
}
