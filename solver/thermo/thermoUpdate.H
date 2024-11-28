Info << "Reading field thermoUpdate\n"<< endl;

forAll(mesh.cells(), i)
{
    // Update the phase fraction due to heat/solute transport
    delta1[i] = alpha10[i] + dFdT[i]*(T[i] - T0[i]);

    delta1[i] = min(max(delta1[i], 0.0), 1.0);

    // Update the phase compositions
    Y2[i] = Y[i] / (1.0 + delta1[i]*(kp.value() - 1.0));

    Y2[i] = min(Y2[i], Yeut.value());
    Y1[i] = kp.value()*Y2[i];
}

delta1.correctBoundaryConditions();
delta2 = scalar(1) - delta1;

Y1.correctBoundaryConditions();
Y2.correctBoundaryConditions();
