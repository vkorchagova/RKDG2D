#include "EdgeInternal.h"


using namespace std;

// ------------------ Constructors & Destructor ----------------



void EdgeInternal::getLocalFluxes(const Flux& flux)
{
    cout << "in internal \n";

    cout << neibCells[0]->number << ' ' << neibCells[1]->number << endl;

    for (int i = 0; i < nGP; ++i)
    {
        numvector<double, 5> solLeft = neibCells[0]->reconstructSolution(gPoints[i]); //alpha coeffs placed in cell
        numvector<double, 5> solRight = neibCells[1]->reconstructSolution(gPoints[i]);


        std::cout << "solLeft" << solLeft << '\n';


        localFluxes[i] = flux.evaluate(solLeft, solRight);
    }
}



