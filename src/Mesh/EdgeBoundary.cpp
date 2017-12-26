#include "EdgeBoundary.h"

using namespace std;

// ------------------ Constructors & Destructor ----------------

void EdgeBoundary::getLocalFluxes(const Flux &flux)
{

    for (int i = 0; i < nGP; ++i)
    {
        numvector<double, 5> solInner = rotate(neibCells[0]->reconstructSolution(gPoints[i]),n);
        numvector<double, 5> solOuter = applyBoundary(solInner);
                        
        localFluxes[i] = flux.evaluate(solInner, solOuter, n);
    }
}
