#include "FluxPh.h"
#include <iostream>

using namespace std;

// ------------------ Constructors & Destructors ----------------


// ------------------ Private class methods --------------------

// ------------------ Public class methods ---------------------



numvector<double, dimPh> FluxPh::evaluate(
        const numvector<double, dimPh>& solLeftRot, const numvector<double, dimPh>& solRightRot,
        const numvector<double, dimGrad>& gradSolLeft, const numvector<double, dimGrad>& gradSolRight
) const
{    
    numvector<double, dimPh> fluxOutward = phs.fluxF(solLeftRot);
    numvector<double, dimPh> fluxInward  = phs.fluxF(solRightRot);

    numvector<double, dimPh> lV = lambdaF_semisum(solLeftRot, solRightRot);
    
    double lambda = max(fabs(lV[0]), fabs(lV[dimPh-1]));

    /*numvector<double, dimPh> f;
    f[0] = fluxOutward[0] <= fluxInward[0] ?  fluxOutward[0] : fluxInward[0];
    f[1] = fluxOutward[1] <= fluxInward[1] ?  fluxOutward[1] : fluxInward[1];
    f[2] = fluxOutward[2] <= fluxInward[2] ?  fluxOutward[2] : fluxInward[2];
    f[3] = fluxOutward[3] <= fluxInward[3] ?  fluxOutward[3] : fluxInward[3];
    f[4] = fluxOutward[4] <= fluxInward[4] ?  fluxOutward[4] : fluxInward[4];
   //fluxOutward) + 0.5 * lambda *  (solLeftRot - solRightRot);*/

    numvector<double, dimPh> f = 0.5 * (fluxInward + fluxOutward) + 0.5 * lambda *  (solLeftRot - solRightRot);
    return f;
}

