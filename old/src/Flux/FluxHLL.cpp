#include "FluxHLL.h"


numvector<double,5> FluxHLL::evaluate( const numvector<double, 5>& solLeft, const numvector<double, 5>& solRight, const Point& n) const
{
    numvector<double, 5> fluxL = problem.fluxF(solLeft);
    numvector<double, 5> fluxR = problem.fluxF(solRight);

    numvector<double, 5> lV = problem.lambdaF(solLeft,solRight);

    if (lV[0] >= 0)
        return inverseRotate(fluxL, n);

    if (lV[4] <= 0)
        return inverseRotate(fluxR, n);

    return inverseRotate((lV[4] * fluxL - lV[0] * fluxR + lV[0] * lV[4] * (solRight - solLeft)) * (1.0 / (lV[4] - lV[0])), n);
}
