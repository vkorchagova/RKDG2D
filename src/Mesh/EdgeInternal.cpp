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
