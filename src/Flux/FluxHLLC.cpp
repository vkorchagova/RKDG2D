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

numvector<double, 5> FluxHLLC::getFStar(const numvector<double, 5>& sol, const numvector<double, 5>& fK, double rhoL, double pK, double SK, double cK, double SStar) const
{
    numvector<double, 5> D = { 0.0, 1.0, 0.0, 0.0, SStar};

    // var 2
    //return ( SStar * (sol*SK - fK) + D * pK ) * (1.0 / (SK - SStar));

    // var 1
    return ( SStar * (sol*SK - fK) + D * SK * (pK + rhoL * cK * (SStar - sol[1] / sol[0])) ) * (1.0 / (SK - SStar));
}

numvector<double,5> FluxHLLC::evaluate( const numvector<double, 5>& solLeft, const numvector<double, 5>& solRight, const Point& n) const
{
    numvector<double, 5> fluxL = problem.fluxF(solLeft);
    numvector<double, 5> fluxR = problem.fluxF(solRight);

    numvector<double, 5> S = problem.lambdaF(solLeft,solRight);

    double SL = min(S[0],S[4]);
    double SR = max(S[0],S[4]);

    if (SL >= 0.0)
        return inverseRotate(fluxL, n);

    if (SR <= 0.0)
        return inverseRotate(fluxR, n);

    double pLeft = problem.getPressure(solLeft);
    double pRight = problem.getPressure(solRight);

    double cLeft = SL - solLeft[1] / solLeft[0];
    double cRight = SR - solRight[1] / solRight[0];

    double sStar = (pRight - pLeft + solLeft[1] * cLeft - solRight[1] * cRight) / \
                   (solLeft[0] * cLeft - solRight[0] * cRight);

//    if (sStar >= 0.0)
//        return inverseRotate(fluxL, n) + SL * inverseRotate(getUStar(solLeft,pLeft,SL,cLeft,sStar) - solLeft, n);

//    return inverseRotate(fluxL, n) + SR * inverseRotate(getUStar(solRight,pRight,SR,cRight,sStar) - solRight, n);

    if (sStar >= 0.0)
        return inverseRotate(getFStar(solLeft,fluxL,solLeft[0],pLeft,SL,cLeft,sStar),n);

    return inverseRotate(getFStar(solRight,fluxR,solLeft[0],pRight,SR,cRight,sStar),n);

    double pLR = 0.5 * (pLeft + pRight + solLeft[0] * cLeft * (sStar - solLeft[1] / solLeft[0]) + \
                                         solRight[0] * cRight * (sStar - solRight[1] / solRight[0]));

    if (sStar >= 0.0)
        return inverseRotate(getFStar(solLeft,fluxL,solLeft[0],pLR,SL,cLeft,sStar),n);

    return inverseRotate(getFStar(solRight,fluxR,solLeft[0],pLR,SR,cRight,sStar),n);


}

