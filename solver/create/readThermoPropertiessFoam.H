IOdictionary thermoProperties
(
    IOobject
    (
        "thermoProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dimensionedScalar Tmelt
(
    "Tmelt",
    dimTemperature,
    thermoProperties.lookup("Tmelt")
);

dimensionedScalar Teut
(
    "Teut",
    dimTemperature,
    thermoProperties.lookup("Teut")
);

dimensionedScalar Tiso
(
    "Tiso",
    dimTemperature,
    thermoProperties.lookup("Tiso")
);

dimensionedScalar Ydelta
(
    "Ydelta",
    dimless,
    thermoProperties.lookup("Ydelta")
);

dimensionedScalar Yeut
(
    "Yeut",
    dimless,
    thermoProperties.lookup("Yeut")
);

dimensionedScalar L
(
    "L",
    dimEnergy/dimMass,
    thermoProperties.lookup("L")
);

// solidus and liquidus slopes
dimensionedScalar msol = (Teut - Tmelt)/Ydelta;
dimensionedScalar mliq = (Teut - Tmelt)/Yeut;

dimensionedScalar kp = Ydelta/Yeut;

const scalar nThermoCorr
(
    pimple.dict().lookupOrDefault<scalar>("nThermoCorr", 100)
);
