# RKDG2D

Runge --- Kutta Discontinuous Galerkin for 2D perfect inviscid gas flows

## Branches

### rectangular

Meshes: rectangular

Riemann solvers: Local Lax --- Fridriechs, HLL, HLLC

Limiters: FinDiff, MUSCL, WENO_S

Indicators: Nowhere, Everywhere, Harten, KXRCF

Boundary conditions: Constant (fixed), Slip, Sine (sine function), Open

Time step: static

Parallelisation: no

### OpenMP

Meshes: unstructured (supports triangular and quadrangle elements), converted from .unv

Riemann solvers: Local Lax --- Fridriechs, HLL, HLLC

Limiters: FinDiff, WENO_S, Riemann WENO_S

Indicators: Nowhere, Everywhere, KXRCF

Boundary conditions: Constant (fixed), Slip, Sine (sine function), Open, Periodic

Time step: static, dynamic (like in OpenFOAM)

Parallelisation: OpenMP

### MPI

Meshes: unstructured (supports triangular and quadrangle elements), converted from .unv

Riemann solvers: Local Lax --- Fridriechs, HLL, HLLC

Limiters: FinDiff, WENO_S, BJ, BJVertex (non-parallel mode)

Indicators: Nowhere, Everywhere, Riemann WENO_S

Boundary conditions: Constant (fixed), Slip, Sine (sine function), Open

Time step: static

Parallelisation: MPI (+OpenMP optionally)
