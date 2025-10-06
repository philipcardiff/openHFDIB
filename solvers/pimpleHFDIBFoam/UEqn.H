MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(U) + fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevReff(U)
 == fvOptions(U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();
fvOptions.constrain(UEqn);

// IB contribution (no 'surface' variable anywhere)
volScalarField rAU(1.0/UEqn.A());
HFDIB.interpolateIB(U, Ui, lambda);
f = lambda * ((1.0/rAU)*Ui - UEqn.H() + fvc::grad(p));

solve(UEqn == -fvc::grad(p) + f);
fvOptions.correct(U);
