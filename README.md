# RKDG2D

Runge --- Kutta Discontinuous Galerkin for 2D perfect inviscid gas flows

Meshes: rectangular
Riemann solvers: Local Lax --- Fridriechs, HLL, HLLC
Limiters: FinDiff, MUSCL, WENO_S
Indicators: Nowhere, Everywhere, Harten, KXRCF
Boundary conditions: Constant (fixed), Slip, Sine (sine function), Open
Time step: static
Parallelisation: no