#include "FluxLLF.h"
#include <iostream>

using namespace std;

// ------------------ Constructors & Destructors ----------------


// ------------------ Private class methods --------------------

// ------------------ Public class methods ---------------------



numvector<double,dimPh> FluxLLF::evaluate( const numvector<double, dimPh>& solLeftRot, const numvector<double, dimPh>& solRightRot) const
{    
    numvector<double, dimPh> fluxOutward = phs.fluxF(solLeftRot);
    numvector<double, dimPh> fluxInward  = phs.fluxF(solRightRot);

    numvector<double, dimPh> lV = lambdaF(solLeftRot, solRightRot);
    
    double lambda = max(fabs(lV[0]), fabs(lV[dimPh-1]));

    return 0.5 * (fluxInward + fluxOutward) + 0.5 * lambda * (solLeftRot - solRightRot);
}

