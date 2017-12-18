#include "EdgeInternal.h"


using namespace std;

// ------------------ Constructors & Destructor ----------------


void EdgeInternal::getLocalFluxes(const Flux& flux)
{
    for (int i = 0; i < nGP; ++i)
    {
        numvector<double, 5> solInner = neibCells[0]->reconstructSolution(gPoints[i]);
        numvector<double, 5> solOuter = neibCells[1]->reconstructSolution(gPoints[i]);

        localFluxes[i] = flux.evaluate(solInner, solOuter, n);
    }
}
