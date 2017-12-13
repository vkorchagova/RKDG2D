#include "EdgeInternal.h"


using namespace std;

// ------------------ Constructors & Destructor ----------------



void EdgeInternal::getLocalFluxes() const
{
    //std::cout << "in internal \n";

    for (int i = 0; i < nGP; ++i)
    {
        numvector<double, 5> solLeft = neibCells[0]->reconstructSolution(gPoints[i]); //alpha coeffs placed in cell
        numvector<double, 5> solRight = neibCells[1]->reconstructSolution(gPoints[i]);


        std::cout << "solLeft" << solLeft << '\n';


        //flux.evaluate(solLeft, solRight);
    }
}



