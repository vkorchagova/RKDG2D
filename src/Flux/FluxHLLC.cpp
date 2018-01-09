#include "FluxHLLC.h"
#include <algorithm>

using namespace std;

numvector<double, 5> FluxHLLC::getUStar (const numvector<double, 5>& sol, double lK, double cK, double lStar) const
{
    double rhoL = sol[0] * cK / (lK - lStar);

    double e = sol[4] / sol[0] + (lStar - sol[1] / sol[0]) * ( lStar + problem.getPressure(sol) / sol[0] / cK);

    numvector<double, 5> iU = { 1.0, lStar, sol[2], sol[3], e};

    return  iU * rhoL;
}

numvector<double,5> FluxHLLC::evaluate( const numvector<double, 5>& solLeft, const numvector<double, 5>& solRight, const Point& n) const
{
    numvector<double, 5> fluxL = inverseRotate(problem.fluxF(solLeft), n);
    numvector<double, 5> fluxR = inverseRotate(problem.fluxF(solRight), n);

    numvector<double, 5> lV = problem.lambdaF(solLeft,solRight);

    if (lV[0] >= 0.0)
        return fluxL;

    if (lV[4] <= 0.0)
        return fluxR;

    double lL = min(lV[0],lV[4]);
    double lR = max(lV[0],lV[4]);

    double pLeft = problem.getPressure(solLeft);
    double pRight = problem.getPressure(solRight);

    double cLeft = lL - solLeft[1] / solLeft[0];
    double cRight = lR - solRight[1] / solRight[0];

    double lStar = (pRight - pLeft + solLeft[1] * cLeft - solRight[1] * cRight) / \
                   (solLeft[0] * cLeft - solRight[0] * cRight);

    if (lStar >= 0.0)
        return fluxL + lL * inverseRotate(getUStar(solLeft,lL,cLeft,lStar) - solLeft, n);

    return fluxR + lR * inverseRotate(getUStar(solRight,lR,cRight,lStar) - solRight, n);
}

