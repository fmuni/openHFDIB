volScalarField rAU(1.0/UEqn.A());
volVectorField HbyA("HbyA", U);
HbyA = rAU*UEqn.H();
surfaceScalarField phiHbyA
(
    "phiHbyA",
    (fvc::interpolate(HbyA) & mesh.Sf())
  + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
);

MRF.makeRelative(phiHbyA);

adjustPhi(phiHbyA, U, p);

surfaceScalarField phiIB(fvc::flux(lambda*rAU*f));

// Non-orthogonal pressure corrector loop
while (piso.correctNonOrthogonal())
{
    // Pressure corrector

    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAU, p) ==  fvc::div(phiHbyA)
                                 
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve();

    if (piso.finalNonOrthogonalIter())
    {
        phi = phiHbyA + phiIB - pEqn.flux();
    }
}


U = HbyA - rAU*fvc::grad(p) + rAU*f;
U.correctBoundaryConditions();
//Update imposed field
rUA = 1.0/rAU;
rUA.correctBoundaryConditions();
HFDIB.interpolateIB(U,Ui,lambda);
f = surface*(rUA*Ui - UEqn.H() + fvc::grad(p));

#include "continuityErrs.H"
fvOptions.correct(U);
