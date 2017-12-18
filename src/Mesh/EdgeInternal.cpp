#include "EdgeInternal.h"


using namespace std;

// ------------------ Constructors & Destructor ----------------



void EdgeInternal::getLocalFluxesHor(const Flux& flux)
{
    for (int i = 0; i < nGP; ++i)
    {
        numvector<double, 5> solLeft = neibCells[0]->reconstructSolution(gPoints[i]); //alpha coeffs placed in cell
        numvector<double, 5> solRight = neibCells[1]->reconstructSolution(gPoints[i]);

        localFluxes[i] = flux.evaluate(solLeft, solRight, n);
    }
}


void EdgeInternal::getLocalFluxesVer(const Flux& flux)
{
    for (int i = 0; i < nGP; ++i)
    {
        numvector<double, 5> solLeft = neibCells[0]->reconstructSolution(gPoints[i]); //alpha coeffs placed in cell
        numvector<double, 5> solRight = neibCells[1]->reconstructSolution(gPoints[i]);
        localFluxes[i] = flux.evaluate(solLeft, solRight, n);
    }
}



