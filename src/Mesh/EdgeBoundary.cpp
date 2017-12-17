#include "EdgeBoundary.h"

using namespace std;

// ------------------ Constructors & Destructor ----------------

void EdgeBoundary::getLocalFluxesHor(const Flux &flux)
{

    for (int i = 0; i < nGP; ++i)
    {
        numvector<double, 5> solLeft = neibCells[0]->reconstructSolution(gPoints[i]); 
        numvector<double, 5> solRight = applyBoundary(solLeft);
                        
        localFluxes[i] = flux.evaluateHor(solLeft, solRight);
    }
}

void EdgeBoundary::getLocalFluxesVer(const Flux &flux)
{

    for (int i = 0; i < nGP; ++i)
    {
        numvector<double, 5> solLeft = neibCells[0]->reconstructSolution(gPoints[i]);
        numvector<double, 5> solRight = applyBoundary(solLeft);

        localFluxes[i] = flux.evaluateVer(solLeft, solRight);
    }
}
