#include "createRDeltaT.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);




// ************************************    From    SolidificationFoam     ****************************//
solidificationSystem fluid(U);

Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (                                       // From    SolidificationFoam     ****************************//
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
// ************************************    From    SolidificationFoam     ****************************//






#include "createPhi.H"

Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());
#include "createFaceFractions.H"

const dimensionedScalar& rho1 = mixture.rho1();	//Use a copy instead of reference. Do not require this as we get this from fluid.rho() of sFoam
//const dimensionedScalar& rho1 = fluid.rho().dimensionedInternalField();
const dimensionedScalar& rho2 = mixture.rho2();         // rho2 of air

// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
	IOobject::AUTO_WRITE
    ),
  //  alpha1*rho1 + alpha2*rho2
   alpha1*fluid.rho() + alpha2*rho2
);
rho.oldTime();

#include "createAuxilaryFieldssFoam.H"    // added from solidificationFoam
#include "createScalarFields.H"

volScalarField mu
(
    IOobject
    (
        "mu",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
    mixture.mu()
);

// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);




// ************************************    From    SolidificationFoam     ****************************//
//#include "compressibleCreatePhi.H"

// Kinematic density for buoyancy force
volScalarField rhok
(
    IOobject                                           // From    SolidificationFoam     ****************************//
    (
        "rhok",
        runTime.timeName(),
        mesh
    ),
    fluid.rho2()
   *(1.0 - fluid.betaT()*(T - fluid.TRef()) - fluid.betaC()*(Y2 - fluid.CRef()))
);
// ************************************    From    SolidificationFoam     ****************************//


//Some changes required below rho ==rhok


// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture)
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

//Check p and p_rgh initialization in incompressible vs compressible.
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

volVectorField Um
(
    IOobject
    (
        "Um",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    U
);

volScalarField alpham
(
    IOobject
    (
        "alpham",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alpha1
);

#include "createMRF.H"

volScalarField surfaceCell
(
    IOobject
    (
        "surfaceCell",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

#include "createCustomDict.H"
//#include "createScalarFields.H"
#include "createScalarSource.H"
#include "createIsoAdvection.H"

#include "createFvOptions.H"
