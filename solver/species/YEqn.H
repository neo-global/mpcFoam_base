
//This replaces CEqn.H. Discuss in case of doubt.

dictionary& customDict = mesh.lookupObjectRef<IOdictionary>("values");
word solveScalar(customDict.lookup("solveScalar"));
if(solveScalar=="true")
{
	//get phase reference
        volScalarField& alpha = mesh.lookupObjectRef<volScalarField>
	(
	    IOobject::groupName("alpha", mixture.alpha1().group())
	);

        surfaceScalarField& alphaf = mesh.lookupObjectRef<surfaceScalarField>
	(
	    IOobject::groupName("alphaf", mixture.alpha1().group())
	);

	//perform phase bounding of convection-diffusion
	Dt = (turbulence->nut()/Sct + DY);
	surfaceScalarField Dtf = alphaf*fvc::interpolate(Dt);

	word useDtf(customDict.lookup("useDtf"));
	if(useDtf=="true")
	{
	    //Dtf *= alphaf;
	    Dtf=Dtf*(alphaf*rho1_B);
	}

        //surfaceScalarField alphaPhi = fvc::flux(U)*alphaf;
        surfaceScalarField alphaPhi = rhoPhi*alphaf;
        

	word interfaceMethod
	(
	    mesh.solverDict("alpha").lookupOrDefault<word>("interfaceMethod", "isoAdvector")
	);	
	if(interfaceMethod=="isoAdvector")
	{
	     alphaPhi = advector.alphaPhi();
	}


        //Scheme
        word divScheme("div(phi,Y)");

	//Solve YEqn
        
        /* // Relative composition for advection/diffusion source terms
              volScalarField CRel("CRel", C2 - C);

         fvScalarMatrix CEqn
        (
        fvm::ddt(rho, C)                              from SolidificationFoam
      + fvm::div(phi, C)
      - fvm::laplacian(alpha2*rho2*D2, C)
    ==
        fvc::laplacian(alpha2*rho2*D2, CRel)
      - fvc::div(phi, CRel)
    );

        */   
        
        volScalarField YRel("YRel", Y2 - Y);
           
        fvScalarMatrix YEqn
        (
             fvm::ddt(rho,Y)
           + fvm::div(alphaPhi, Y, divScheme) //replace this with alphaPhi1Un
           - fvm::laplacian(Dtf, Y)
           ==
           fvc::laplacian(Dtf, YRel)
           - fvc::div(alphaPhi, YRel)
           //  ==
           // fvOptions(rho, Y)
           // - fvm::Sp(S_Y, Y)
           // + S_Y*Ystar*Yg
        );

        YEqn.relax();
        fvOptions.constrain(YEqn);
        YEqn.solve(mesh.solver("Y"));

        fvOptions.correct(Y);

        //write fields for debugging
	//Comment to save space
        if(mesh.time().writeTime())
        {
            volVectorField alphaPhiVol
            (
        	IOobject
        	(
                    "yEqnFlux",
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
        	),
        	fvc::reconstruct(alphaPhi)
            );

            alphaPhiVol.write();
            Dt.write();
            Y.write();
        }
}
