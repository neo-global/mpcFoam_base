dictionary& customDict = mesh.lookupObjectRef<IOdictionary>("values");
word solveScalar(customDict.lookup("solveScalar"));
if(solveScalar=="true")
{
	//get phase reference
        volScalarField& alpha = mesh.lookupObjectRef<volScalarField>(IOobject::groupName("alpha",mixture.alpha1().group()));
        surfaceScalarField& alphaf = mesh.lookupObjectRef<surfaceScalarField>(IOobject::groupName("alphaf",mixture.alpha1().group()));


	//perform phase bounding of convection-diffusion
        surfaceScalarField alphaPhi = advector.alphaPhi(); //check what combination works between U,rho,alpha. Take hint from rhoPhi.
	
	word useDtf(customDict.lookup("useDtf"));

        scalar Sct (readScalar(customDict.lookup("Sct")));
	dimensionedScalar DY
	(
	    "DY",
	    mu.dimensions()/rho.dimensions(),
	    readScalar
	    (
		customDict.lookup("DY")
	    )
	);

        Dt = (turbulence->nut()/Sct+DY);
	surfaceScalarField Dtf = fvc::interpolate(Dt);

	if(useDtf=="true")
	{
	    Dtf *= alphaf;
	}

	//Create sources

	//interface source
    if
    (
        mesh.foundObject<volVectorField>(IOobject::groupName("interfaceNormal",mixture.alpha1().group()))
    )
    {
        volVectorField& area = mesh.lookupObjectRef<volVectorField>(IOobject::groupName("interfaceNormal",mixture.alpha1().group()));

        S_Y.primitiveFieldRef() = 0.4
                                *surfaceSource
                                *Foam::pow(DY, 0.5)
                                *Foam::pow(mixture.nu(),-0.25)
                                *Foam::pow(turbulence->epsilon(),0.25)
                                *mag(area)/mesh.V();
    }

    //dispersed source
/*
    YbubSource = dimensionedScalar("zero", dimless/dimTime, 0.0);
    if(bubbleSource)
    {
        kinematicCloud.massCollection(Yexp, YbubSource);
    }
*/
    //gas-side-reduced-fraction
    scalar YgStar (readScalar(customDict.lookup("YgStar", 1.0)));

    //Scheme
    word divScheme("div(phi,Y)");

    //Solve explicit
/*
 * fvScalarMatrix YexpEqn
    (
         fvm::ddt(Yexp)
       + fvm::div(alphaPhi, Yexp, divScheme) //replace this with alphaPhi1Un
       - fvm::laplacian(Dtf, Yexp)
         ==
        //fvOptions(rho, Y)
        S_Y*(Ystar*YgStar-Yexp)
      + YbubSource
    );

    YexpEqn.relax();
    fvOptions.constrain(YexpEqn);

    YexpEqn.solve(mesh.solver("Y"));
    fvOptions.correct(Yexp);

*/

//basic implicit addition to coefficient for Yimplicit.
/*    if(bubbleSource)
	{
        S_Y.primitiveFieldRef() += 0.4
                                  *Foam::pow(DY,0.5)
                                  *Foam::pow(mixture.nu(),-0.25)
                                  *Foam::pow(turbulence->epsilon(),0.25)
                                  *6./(kinematicCloud.Dmax()+VSMALL)
                                  *kinematicCloud.theta();
	}
*/
	//Solve YEqn

        fvScalarMatrix YEqn
        (
             fvm::ddt(Y)
           + fvm::div(alphaPhi, Y, divScheme) //replace this with alphaPhi1Un
           - fvm::laplacian(Dtf, Y)
             ==
             //fvOptions(rho, Y)
           - fvm::Sp(S_Y, Y)
           + S_Y*Ystar*YgStar
        );

        YEqn.relax();
        fvOptions.constrain(YEqn);
        YEqn.solve(mesh.solver("Y"));
        fvOptions.correct(Y);


	//write fields for debugging
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
        	fvc::reconstruct(alphaPhi1Un)
        );

        alphaPhiVol.write();
        Dt.write();
        Y.write();
 //       Yexp.write();
  //      YbubSource.write();
    }
}
