#include "EdgeBoundaryInfty.h"

void EdgeBoundaryInfty::setLocalFluxes()
{
    for (int i = 0; i < nGP; ++i)
    {
        numvector<double, 5> solLeft = neibCells[0].reconstructSolution(gPoints[i]); //alpha coeffs placed in cell
        numvector<double, 5> solRight = applyBoundary(solLeft);

        flux.evaluate(solLeft, solRight);
    }
}

numvector<double, 5> EdgeBoundaryInfty::applyBoundary()
{
    return infty;
}
