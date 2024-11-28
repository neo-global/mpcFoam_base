/*---------------------------------------------------------------------------*\
              Original work | Copyright (C) 2016-2017 DHI
              Modified work | Copyright (C) 2016-2017 OpenCFD Ltd.
              Modified work | Copyright (C) 2017-2018 Johan Roenby
-------------------------------------------------------------------------------

License
    This file is part of IsoAdvector, which is an unofficial extension to
    OpenFOAM.

    IsoAdvector is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    IsoAdvector is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with IsoAdvector. If not, see <http://www.gnu.org/licenses/>.

Application
    interFlow

Description
    Solver for 2 incompressible, isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach,
    with optional mesh motion and mesh topology changes including adaptive
    re-meshing.

\*---------------------------------------------------------------------------*/

#include "isoAdvection.H"
#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
//#include "fieldToCell.H"
//#include "cellSet.H"

//***************************************from SolidificationFoam***************************//
#include "create/solidificationSystem/solidificationSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "create/createFields.H"
    #include "create/createFieldRefssFoam.H"        // added from SolidificationFoam
   // #include "createTimeControls.H"     // added from SolidificationFoam Not Sure about this it is required or not.
    #include "create/readThermoPropertiessFoam.H"  // added from SolidificationFoam
    #include "createFunctionObjects.H"
//    #include "createTurbulenceWakeSource.H"
    #include "createAlphaFluxes.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"
    
    
    
   // #include "setRootCase.H"
   // #include "createTime.H"
   // #include "createMesh.H"
    //#include "createControl.H"
   // #include "createFields.H"
   //#include "createFieldRefs.H"
  //   #include "createTimeControls.H"
    //#include "readThermoProperties.H"

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

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {


	    if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
		    mesh.update();

       	        if (mesh.changing())
                {
    	   	    // Do not apply previous time-step mesh compression flux
                    // if the mesh topology changed
                    if (mesh.topoChanging())
                    {
                        talphaPhi1Corr0.clear();
                    }

                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;

                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);

                        mixture.correct();
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }
            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"

            mixture.correct();

            #include "thermoLoopsFoam.H" // Added
            #include "UEqn.H"
            
            

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }

           // #include "YEqn.H"

        }

	runTime.write();

      	Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
   // #include "dimensionedType.H"

// Example field, e.g., rho or phi
/*Info << "Dimensions of rho: " << rho.dimensions() << endl;
Info << "Dimensions of phi: " << phi.dimensions() << endl;*/

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
