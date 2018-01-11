#include "FluxLLF.h"
#include <iostream>

using namespace std;

// ------------------ Constructors & Destructors ----------------


// ------------------ Private class methods --------------------

// ------------------ Public class methods ---------------------



numvector<double,5> FluxLLF::evaluate( const numvector<double, 5>& solLeft, const numvector<double, 5>& solRight, const Point& n) const
{    
    numvector<double, 5> fluxOutward = problem.fluxF(solLeft);
    numvector<double, 5> fluxInward = problem.fluxF(solRight);

    numvector<double, 5> lV = problem.lambdaF(solLeft,solRight);
    
    double lambda = max(fabs(lV[0]), fabs(lV[4]));

    return inverseRotate(0.5 * (fluxInward + fluxOutward) + 0.5 * lambda * (solLeft - solRight), n);
}
