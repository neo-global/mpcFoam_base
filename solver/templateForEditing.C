/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    mpcFoam

Description
    mpc stands for multiphase phase change. The solver attempts to merge 
    multiphase and sharp-interface capturing functionality of VoF with phase
    change dynamics of alloys including phase transformation, seggregation of
    chemical species, and associated energy dynamics. Model is compatible with
    turbulence modeling as well as choice of VoF strategy. The latter 
    is incorporated through runTime selectable VoF strategy and relies primarily
    on the isoAdvector geometric VoF method from the TwoPhaseFlow library adapted 
    from Roenby et al.
    
    Summary of modules:
    - Solidification and segregation within dense phase: 
      OpenFOAM-solidification by XYZ according to Inconera et al
    - Geometric VoF:
      TwoPhaseFlow by DLR-RY according to Schuler et al
    - Species transport in multiphase domain:
      Mishra et al
    - Base solver:
      Compressible/Incompressible twoPhaseTransport.
   
      Copyright: 
	Rishikesh, Mrigank, Amrita, Biswajit and Amarendra Kumar Singh,
	IIT Kanpur 

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "fvcSmooth.H"
#include "CorrectPhi.H"
#include "dynamicFvMesh.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "pimpleControl.H"
#include "subCycle.H"
//compressible/incompressible
#include "turbulentTransportModel.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
//VoF
#include "CMULES.H"
#include "isoAdvection.H"
//solidification


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
//Check which one works
//    #include "setRootCase.H"
//    #include "setRootCaseLists.H"

    //Time and Space (Mesh)
    #include "createTime.H"
    #include "createRDeltaT.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"

    //createFields
    #include "createFields.H"
    #include "createAlphaFluxes.H" //Used by speciesTransport within VoF framework

    //Flux calculation
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"

    turbulence->validate();

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "alphaCourantNo.H"
            #include "setDeltaT.H"
        }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;
	{ 
	    //This section needs to be adapted
	    //according to compressible VoF and solidification
        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
	}

	//Routine checks...
	//alphaEqn
	//thermoLoop

	// This part is critical to synchronize 
	// VoF, species and Solidificaiton
	{
	    //rhoDense = get from solidification
	    //muDense = get from solidification
	    //rho = alpha*rhoDense + (1.0-alpha)*rho2
	    //mu = alpha*muDense + (1.0-alpha)*mu2
	}
	//UEqn
	//pEqn
	//turbulence
        }

    }

    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
