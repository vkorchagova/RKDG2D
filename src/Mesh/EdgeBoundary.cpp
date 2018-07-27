#include "EdgeBoundary.h"

using namespace std;

// ------------------ Constructors & Destructor ----------------

// ------------------ Public class methods ---------------------

void EdgeBoundary::getLocalFluxes(const Flux &flux)
{

    for (int i = 0; i < nGP; ++i)
    {
        numvector<double, 5> solInner = rotate(neibCells[0]->reconstructSolution(gPoints[i]),n);
        numvector<double, 5> solOuter = bc->applyBoundary(solInner,n);
                        
        localFluxes[i] = flux.evaluate(solInner, solOuter, n);
    }
}

void EdgeBoundary::getMaxUL()
{
    double uMax = 0.0;

    for (size_t i = 0; i < neibCells.size(); ++i)
    {
        numvector<double, 5> solLeft = rotate(neibCells[0]->reconstructSolution(nodes[i]),n);
        numvector<double, 5> solRight = bc->applyBoundary(solLeft);

        numvector<double, 5> lambda = neibCells[0]->problem.lambdaF(solLeft,solRight);

        uMax = max( max(fabs(lambda[0]), fabs(lambda[4])), uMax);
    }

    uMaxL = uMax * getLength();
}
