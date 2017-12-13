#include "EdgeBoundary.h"

using namespace std;

// ------------------ Constructors & Destructor ----------------

void EdgeBoundary::getLocalFluxes(const Flux &flux)
{
    cout << "in boundary edge\n";

    for (int i = 0; i < nGP; ++i)
    {
        numvector<double, 5> solLeft = neibCells[0]->reconstructSolution(gPoints[i]); 
        numvector<double, 5> solRight = applyBoundary(solLeft);
                        
        localFluxes[i] = flux.evaluate(solLeft, solRight);
    }
}
