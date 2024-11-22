# mpcFoam_base
multiphase phase change solver, incompressible

Compile in the following order:
1. interpolationSchemes - these are some generic libraries for general compatibility
2. twoPhaseFlow/src/VoF - this is the specialized VOF code required for better species transport ( CEqn/YEqn )
3. solver - the main solver, currently has two phase incompressible VOF with species transport (YEqn)
