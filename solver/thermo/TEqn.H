{
    // Update mixture properties
    volScalarField Cp(fluid.Cp());
    volScalarField kappa("kappa", fluid.kappa());
   
    // Semi-implicit handling of eutectic reaction
    volScalarField Sp = rho1_B*L*dFdT/runTime.deltaT();  //[1,-1,-3,1,0,0,0]
    volScalarField Sc = -Sp*T0 - fvc::ddt(rho1_B*L, delta2);//    [1,-1,-3,0,0,0,0]

    fvScalarMatrix TEqn
    (
        fvm::ddt(rho*Cp, T)
      + fvm::div(rhoPhi*fvc::interpolate(Cp), T, "div(phi,T)") // dimension[1,1,-3,0,0,0]
      - fvm::laplacian(kappa, T)
     ==
        fvm::Sp(Sp, T) + Sc
      - L*fvc::div(rhoPhi)               //LatentHeatTerm    // dimension[1,1,-3,0,0,0]
    );

    TEqn.relax();
    TEqn.solve();
}


