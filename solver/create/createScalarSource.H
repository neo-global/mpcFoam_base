scalar Yst (readScalar(customDict.lookup("Ystar")));
scalar surfaceSource(readScalar(customDict.lookup("surfaceSource01")));
scalar bubbleSource(readScalar(customDict.lookup("bubbleSource01")));
scalar convectionBound(readScalar(customDict.lookup("convectionBound")));
scalar diffusionBound(readScalar(customDict.lookup("diffusionBound")));
dimensionedScalar DY
(
    "DY",
    mixture.nu()().dimensions(),
    readScalar(customDict.lookup("DY"))
);

/*
volScalarField Vcell
(
    IOobject
    (
	"Vcell",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimVolume, 1.0)
);
Vcell.primitiveFieldRef() = mesh.V();
*/


volScalarField S_Y
(
    IOobject
    (
	"S_Y",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0,0,-1,0,0,0,0), 0.0)
);



volScalarField S_Y1
(
    IOobject
    (
	"S_Y2",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
    S_Y
);

volScalarField S_Y2
(
    IOobject
    (
	"S_Y2",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
    S_Y
);

volScalarField Ystar
(
    IOobject
    (
	"Ystar",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Ystar", dimless, Yst)
);
