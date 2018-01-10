#include "FluxHLLC.h"
#include <algorithm>

using namespace std;

numvector<double, 5> FluxHLLC::getUStar (const numvector<double, 5>& sol, double pK, double SK, double cK, double SStar) const
{
    double mult = sol[0] * cK / (SK - SStar);

    double e = sol[4] / sol[0] + (SStar - sol[1] / sol[0]) * ( SStar + pK / sol[0] / cK);

    numvector<double, 5> iU = { 1.0, SStar, sol[2], sol[3], e};

    return  iU * mult;
}

numvector<double,5> FluxHLLC::evaluate( const numvector<double, 5>& solLeft, const numvector<double, 5>& solRight, const Point& n) const
{
    numvector<double, 5> fluxL = inverseRotate(problem.fluxF(solLeft), n);
    numvector<double, 5> fluxR = inverseRotate(problem.fluxF(solRight), n);

    numvector<double, 5> S = problem.lambdaF(solLeft,solRight);

    double SL = min(S[0],S[4]);
    double SR = max(S[0],S[4]);

    if (SL >= 0.0)
        return fluxL;

    if (SR <= 0.0)
        return fluxR;

    double pLeft = problem.getPressure(solLeft);
    double pRight = problem.getPressure(solRight);

    double cLeft = SL - solLeft[1] / solLeft[0];
    double cRight = SR - solRight[1] / solRight[0];

    double sStar = (pRight - pLeft + solLeft[1] * cLeft - solRight[1] * cRight) / \
                   (solLeft[0] * cLeft - solRight[0] * cRight);

    if (sStar >= 0.0)
        return fluxL + SL * inverseRotate(getUStar(solLeft,pLeft,SL,cLeft,sStar) - solLeft, n);

    return fluxR + SR * inverseRotate(getUStar(solRight,pRight,SR,cRight,sStar) - solRight, n);
}

