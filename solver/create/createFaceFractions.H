// Used by the VoF library to create accurate face-level data of 
// phases and interface.

surfaceScalarField alpha1f
(
    IOobject
    (
	IOobject::groupName("alphaf", mixture.alpha1().group()),
	runTime.timeName(),
	mesh,
	IOobject::READ_IF_PRESENT,
	IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);
/*
surfaceScalarField alpha2f
(
    IOobject
    (
	IOobject::groupName("alphaf", mixture.alpha1().group()),
	runTime.timeName(),
	mesh,
	IOobject::READ_IF_PRESENT,
	IOobject::AUTO_WRITE
    ),
    scalar(1) - alpha1f
);
*/
