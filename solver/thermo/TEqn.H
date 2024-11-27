{
    // Update mixture properties
    volScalarField Cp(fluid.Cp());
    volScalarField kappa("kappa", fluid.kappa());
   
    // Semi-implicit handling of eutectic reaction
    volScalarField Sp = rho1_B*L*dFdT/runTime.deltaT();
    volScalarField Sc = -Sp*T0 - fvc::ddt(rho1_B*L, delta2);

    fvScalarMatrix TEqn
    (
        fvm::ddt(rho*Cp, T)
      + fvm::div(rhoPhi*fvc::interpolate(Cp), T, "div(phi,T)")
      - fvm::laplacian(kappa, T)
     ==
        fvm::Sp(Sp, T) + Sc
      - L*fvc::div(rhoPhi)               //LatentHeatTerm
    );

    TEqn.relax();
    TEqn.solve();
}


