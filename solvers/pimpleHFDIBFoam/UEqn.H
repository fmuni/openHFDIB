// Solve the Momentum equation
forAll(Ui,cellI)
{
    Ui[cellI].component(0)=0.0;
    if( lambda[cellI]>0.001 ) surface[cellI] =1.0;
    else surface[cellI] =0.0;

}

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
volScalarField rAU(1.0/UEqn.A());

HFDIB.interpolateIB(U,Ui,lambda);
f = surface*((1.0/rAU)*Ui - UEqn.H() + fvc::grad(p) );

solve(UEqn == -fvc::grad(p) + f);

fvOptions.correct(U);

