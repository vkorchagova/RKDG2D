#include "FluxLLF.h"
#include <iostream>

using namespace std;

// ------------------ Constructors & Destructors ----------------


// ------------------ Private class methods --------------------

// ------------------ Public class methods ---------------------



numvector<double,5> FluxLLF::evaluate( const numvector<double, 5>& solLeft, const numvector<double, 5>& solRight, const Point& n) const
{

    //numvector<double, 5> fluxOutward = n.x() * problem.fluxF(solInner) + n.y() * problem.fluxG(solInner);
    //numvector<double, 5> fluxInward = n.x() * problem.fluxF(solOuter) + n.y() * problem.fluxG(solOuter);
    
    numvector<double, 5> fluxOutward = inverseRotate(problem.fluxF(solLeft), n);
    numvector<double, 5> fluxInward = inverseRotate(problem.fluxF(solRight), n);

    numvector<double, 5> lV = problem.lambdaF(solLeft,solRight);
    
    double lambda = max(fabs(lV[0]), fabs(lV[4]));

    return 0.5 * (fluxInward + fluxOutward) + 0.5 * lambda * (solLeft - solRight);
}
