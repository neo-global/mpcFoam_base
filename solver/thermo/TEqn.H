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


/*Info<< "Reading field Div\n" << endl;
volScalarField Div =fvc::div(rhoPhi*fvc::interpolate(Cp), T, "div(phi,T)");

//dimensionedScalar Div =((phi*fvc::interpolate(Cp))*T);
surfaceScalarField INCP=fvc::interpolate(Cp);*/

    
    
    
    
    
    
    
//#include "dimensionedType.H"
//Info << "Reading field p_rgh\n" << rho.dimensions()<< endl;
//Info << "Reading field kappa\n" << kappa.dimensions()<< endl;
/*Info << "Reading field Sp\n" << Sp.dimensions()<< endl;
Info << "Reading field Sc\n" << Sc.dimensions()<< endl;
Info << "Reading field rhoPhi\n" << rhoPhi.dimensions()<< endl;
Info << "Reading field L\n" << L.dimensions()<< endl;
//Info << "Reading field L1\n" << L1.dimensions()<< endl;
Info << "Reading field Div\n" << Div.dimensions()<< endl;
Info << "Reading field Cp\n" << Cp.dimensions()<< endl;
Info << "Reading field INCP\n" << INCP.dimensions()<< endl;*/
    TEqn.relax();
    TEqn.solve();
}


