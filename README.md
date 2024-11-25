# mpcFoam_base
multiphase phase change solver, incompressible

Compile in the following order:

0. Apply triSurface patch to $FOAM_SRC/triSurface/triSurface.H
```
cp $FOAM_SRC/triSurface/triSurface.H $FOAM_SRC/triSurface/triSurface.backup
patch -u $FOAM_SRC/triSurface/triSurface.H -i triSurface/triSurfaceMod.patch 
```
1. compile interpolationSchemes - these are some generic libraries for general compatibility
2. compile twoPhaseFlow/src/VoF - this is the specialized VOF code required for better species transport ( CEqn/YEqn )
3. compile solver - the main solver, currently has two phase incompressible VOF with species transport (YEqn)
