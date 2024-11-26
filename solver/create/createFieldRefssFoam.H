volScalarField& delta1(fluid.alpha1());
volScalarField& delta2(fluid.alpha2());

const dimensionedScalar& rho1_B = fluid.rho2();

const dimensionedScalar& mu2 = fluid.mu2();

const dimensionedScalar& D2 = fluid.D2();

const dimensionedScalar& DAS = fluid.DAS();         // what is this?

const dimensionedScalar& betaT = fluid.betaT();
const dimensionedScalar& betaC = fluid.betaC();

const dimensionedScalar& TRef = fluid.TRef();
const dimensionedScalar& CRef = fluid.CRef();
