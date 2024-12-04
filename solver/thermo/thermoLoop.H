for (int tCorr=0; tCorr<nThermoCorr; tCorr++)
{
    // Store previous values for linearized update
    volScalarField alpha10(delta1);
    volScalarField T0(T);

    // Update the freezing range
    volScalarField Tliq(min(max(Tmelt + mliq*Y, Teut), Tmelt));
    volScalarField Tsol(min(max(Tmelt + msol*Y, Teut), Tmelt));

    #include "thermo/thermoSource.H"
    #include "thermo/TEqn.H"
    #include "thermo/thermoUpdate.H"
    #include "thermo/YEqn.H"

    // Check for phase fraction convergence
    scalarField residualField
    (
        mag(delta1.internalField() - alpha10.internalField())
    );

    scalar residual(gMax(residualField));

    Info<< "THERMO: iteration " << (tCorr + 1)
        << " residual: " << residual  << endl;

    if (residual < 1e-4)
    {
        break;
    }



}
