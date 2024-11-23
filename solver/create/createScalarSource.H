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

