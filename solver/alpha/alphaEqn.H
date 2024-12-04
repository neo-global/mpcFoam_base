// Looking up which alpha1 solver to use
const dictionary& alpha1Dict(mesh.solverDict(alpha1.name()));

word interfaceMethod
(
    alpha1Dict.lookupOrDefault<word>("interfaceMethod", "isoAdvector")
);

// Advancing alpha1 field
if (interfaceMethod == "MULES")
{
    #include "alphaEqnMULES.H"
    advector.construct();
}
else if (interfaceMethod == "isoAdvector")
{
    // Updating alpha1
    advector.advect();
    #include "rhofs.H"
    rhoPhi = advector.getRhoPhi(rho1f, rho2f);

    alpha2 = 1.0 - alpha1;
    mixture.correct();

    Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  1 - Max(" << alpha1.name() << ") = " << 1 - max(alpha1).value()
        << endl;

    advector.construct();
}
