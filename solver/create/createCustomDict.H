IOdictionary customDict
(
    IOobject
    (
	"values",
	runTime.constant(),
	mesh,
	IOobject::MUST_READ,
	IOobject::NO_WRITE
    )
);

dimensionedScalar Sct
(
    "Sct",
    dimless,
    readScalar(customDict.lookup("Sct"))
);


