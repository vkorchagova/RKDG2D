#include "FluxLLF.h"
#include <iostream>

using namespace std;

// ------------------ Constructors & Destructors ----------------


// ------------------ Private class methods --------------------

// ------------------ Public class methods ---------------------

numvector<double,5> FluxLLF::evaluate( const numvector<double, 5>& solInner, const numvector<double, 5>& solOuter, const Point& n) const
{

    numvector<double,5> fluxOutward = n.x() * problem->fluxF(solInner) + n.y() * problem->fluxG(solInner);
    numvector<double,5> fluxInward = n.x() * problem->fluxF(solOuter) + n.y() * problem->fluxG(solOuter);

    double lambda = n.x() * problem->lambdaF(solInner,solOuter)[4] + n.y() * problem->lambdaG(solInner,solOuter)[4];

    return 0.5 * (fluxInward + fluxOutward) + 0.5 * lambda * (solInner - solOuter);
}
