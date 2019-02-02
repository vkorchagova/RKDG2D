#include "FluxHLL.h"


numvector<double,dimPh> FluxHLL::evaluate( 
    const numvector<double, dimPh>& solLeft, 
    const numvector<double, dimPh>& solRight
) const
{
    numvector<double, dimPh> fluxL = phs.fluxF(solLeft);
    numvector<double, dimPh> fluxR = phs.fluxF(solRight);

    numvector<double, dimPh> lV = lambdaF(solLeft,solRight);

    if (lV[0] >= 0)
        return fluxL;

    if (lV[4] <= 0)
        return fluxR;

    return (lV[4] * fluxL - lV[0] * fluxR + lV[0] * lV[4] * (solRight - solLeft)) * (1.0 / (lV[4] - lV[0]));
}
