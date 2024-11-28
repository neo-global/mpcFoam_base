
//This replaces CEqn.H. Discuss in case of doubt.
Info << "Reading field YEqn\n"<< endl;
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
	    Dtf=Dtf*(alphaf);
	}

        //surfaceScalarField alphaPhi = fvc::flux(U)*alphaf;
        surfaceScalarField alphaPhi= rhoPhi*alphaf;
        

	word interfaceMethod
	(
	    mesh.solverDict("alpha").lookupOrDefault<word>("interfaceMethod", "isoAdvector")
	);	
	if(interfaceMethod=="isoAdvector")
	{     
	     Info << "Reading field Advector\n"<< endl;
	     alphaPhi = rho1_B*advector.alphaPhi();
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
     /*  Info<< "Reading field AlphaDiv\n" << endl;
volScalarField AlphaDiv =fvc::div(alphaPhi , Y, divScheme);

Info<< "Reading field Ddt\n" << endl;
volScalarField Ddt =fvc::ddt(rho,Y);

Info<< "Reading field LapL\n" << endl;
volScalarField LapL =fvc::laplacian(Dtf, Y);

Info<< "Reading field LapR\n" << endl;
volScalarField LapR =fvc::laplacian(Dtf, YRel);

Info<< "Reading field DivR\n" << endl;
volScalarField DivR =fvc::div(alphaPhi, YRel);
 

//dimensionedScalar Div =((phi*fvc::interpolate(Cp))*T);
//surfaceScalarField INCP=fvc::interpolate(Cp);*/

    
    
    
    
    
    
    
/*#include "dimensionedType.H"
Info << "Reading field alphaPhi\n" << alphaPhi.dimensions()<< endl;
Info << "Reading field Y\n" << Y.dimensions()<< endl;
Info << "Reading field Dtf\n" << Dtf.dimensions()<< endl;
Info << "Reading field rho1_B\n" << rho1_B.dimensions()<< endl;
Info << "Reading field rhoPhi\n" << rhoPhi.dimensions()<< endl;
//Info << "Reading field L\n" << L.dimensions()<< endl;
//Info << "Reading field L1\n" << L1.dimensions()<< endl;
Info << "Reading field AlphaDiv\n" << AlphaDiv.dimensions()<< endl;
Info << "Reading field Ddt\n" << Ddt.dimensions()<< endl;
Info << "Reading field LapL\n" << LapL.dimensions()<< endl;
Info << "Reading field LapR\n" << LapR.dimensions()<< endl;
Info << "Reading field DivR\n" << DivR.dimensions()<< endl;
//Info << "Reading field INCP\n" << INCP.dimensions()<< endl;*/ 
        
        
        
      
        
           
        fvScalarMatrix YEqn
        (
             fvm::ddt(rho,Y)
           + fvm::div(alphaPhi, Y, divScheme) //replace this with alphaPhi1Un
           - fvm::laplacian(Dtf*rho1_B, Y)
           ==
           fvc::laplacian(Dtf*rho1_B, YRel)
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
        
        
    if(mesh.time().writeTime())
        {
            volScalarField Cp
            (
        	IOobject
        	(
                    "Cp",
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
        	),
        	fluid.Cp()
            );

            Cp.write();
            Dt.write();
            Y.write();
        }
}
