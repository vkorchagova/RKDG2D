#include "EdgeInternal.h"


using namespace std;

// ------------------ Constructors & Destructor ----------------


void EdgeInternal::getLocalFluxes(const Flux& flux)
{
    for (int i = 0; i < nGP; ++i)
    {
        numvector<double, 5> solInner = rotate(neibCells[0]->reconstructSolution(gPoints[i]),n);
        numvector<double, 5> solOuter = rotate(neibCells[1]->reconstructSolution(gPoints[i]),n);

        localFluxes[i] = flux.evaluate(solInner, solOuter, n);
    }
}

void EdgeInternal::getMaxUL()
{
    double uMax = 0.0;

    for (int i = 0; i < 2; ++i)
    {
        numvector<double, 5> solLeft  = rotate(neibCells[0]->reconstructSolution(nodes[i]),n);
        numvector<double, 5> solRight = rotate(neibCells[1]->reconstructSolution(nodes[i]),n);

        numvector<double, 5> lambda = neibCells[0]->problem.lambdaF(solLeft,solRight);

        uMax = max( max(fabs(lambda[0]), fabs(lambda[4])), uMax);
    }

    uMaxL = uMax * getLength();
}
