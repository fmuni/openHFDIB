// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

fvVectorMatrix UEqn
(
    fvm::ddt(U) + fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevReff(U)
 ==
    fvOptions(U)
);

UEqn.relax();

fvOptions.constrain(UEqn);

//if (piso.momentumPredictor())
//{
    volScalarField rUA(UEqn.A());
    rUA.correctBoundaryConditions();
    HFDIB.interpolateIB(U,Ui,lambda);
    f = surface*(rUA*Ui - UEqn.H() +fvc::grad(p) );
    
    solve(UEqn == -fvc::grad(p) + f );
    
    fvOptions.correct(U);
//}
