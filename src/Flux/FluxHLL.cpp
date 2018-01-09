#include "FluxHLL.h"


numvector<double,5> FluxHLL::evaluate( const numvector<double, 5>& solLeft, const numvector<double, 5>& solRight, const Point& n) const
{
    numvector<double, 5> fluxL = inverseRotate(problem.fluxF(solLeft), n);
    numvector<double, 5> fluxR = inverseRotate(problem.fluxF(solRight), n);

    numvector<double, 5> lV = problem.lambdaF(solLeft,solRight);

    if (lV[0] >= 0)
        return fluxL;

    if (lV[4] <= 0)
        return fluxR;

    return (lV[4] * fluxL - lV[0] * fluxR + lV[0] * lV[4] * inverseRotate(solRight - solLeft, n)) * (1.0 / (lV[4] - lV[0]));
}
