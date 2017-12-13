#include "EdgeBoundary.h"

using namespace std;

// ------------------ Constructors & Destructor ----------------

void EdgeBoundary::getLocalFluxes() const
{
    for (int i = 0; i < nGP; ++i)
    {
        numvector<double, 5> solLeft = neibCells[0]->reconstructSolution(gPoints[i]); 
        numvector<double, 5> solRight = applyBoundary(solLeft);
                        
        //flux.evaluate(solLeft, solRight);
    }
}
