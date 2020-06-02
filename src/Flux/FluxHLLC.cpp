#include "FluxHLLC.h"
#include <algorithm>

using namespace std;

numvector<double, dimPh> FluxHLLC::getUStar (
    const numvector<double, dimPh>& sol, 
    double pK, 
    double SK, 
    double cK, 
    double SStar
) const
{
    double mult = sol[0] * cK / (SK - SStar);

    double e = sol[4] / sol[0] + (SStar - sol[1] / sol[0]) * ( SStar + pK / sol[0] / cK);

    numvector<double, dimPh> iU = { 1.0, SStar, sol[2], sol[3], e};

    return  iU * mult;
}

numvector<double, dimPh> FluxHLLC::getFStar(
    const numvector<double, dimPh>& sol, 
    const numvector<double, dimPh>& fK, 
    double rhoL, 
    double pK, 
    double SK, 
    double cK, 
    double SStar
) const
{
    numvector<double, dimPh> D = { 0.0, 1.0, 0.0, 0.0, SStar};

    // var 2
//    return ( SStar * (sol*SK - fK) + D * pK * SK ) * (1.0 / (SK - SStar));

    // var 1
    return ( SStar * (sol*SK - fK) + D * SK * (pK + rhoL * cK * (SStar - sol[1] / sol[0])) ) * (1.0 / (SK - SStar));
}

numvector<double,dimPh> FluxHLLC::evaluate( 
    const numvector<double, dimPh>& solLeft, 
    const numvector<double, dimPh>& solRight 
) const
{
    numvector<double, dimPh> fluxL = phs.fluxF(solLeft);
    numvector<double, dimPh> fluxR = phs.fluxF(solRight);

//    if (fabs (fluxL[2]) > 1e-6 || fabs (fluxR[2]) > 1e-6)
//        cout << fluxL[2] << ' ' << fluxR[2] << endl;

    numvector<double, dimPh> S = lambdaF(solLeft,solRight);

    double SL = min(S[0],S[4]);
    double SR = max(S[0],S[4]);

    numvector<double, dimPh> res(0.0);

    if (SL > 0.0)
        return fluxL;

    if (SR < 0.0)
        return fluxR;

    double pLeft = phs.getPressure(solLeft);
    double pRight = phs.getPressure(solRight);

    double cLeft = SL - solLeft[1] / solLeft[0];
    double cRight = SR - solRight[1] / solRight[0];

    double sStar = (pRight - pLeft + solLeft[1] * cLeft - solRight[1] * cRight) / \
                   (solLeft[0] * cLeft - solRight[0] * cRight);

    if (sStar > 0.0)
        res = fluxL + (getUStar(solLeft,pLeft,SL,cLeft,sStar) - solLeft) * SL;
    else
        res = fluxR + (getUStar(solRight,pRight,SR,cRight,sStar) - solRight) * SR;

//    return res;

    // var 1

    if (sStar > 0.0)
        return getFStar(solLeft,fluxL,solLeft[0],pLeft,SL,cLeft,sStar);

    return getFStar(solRight,fluxR,solLeft[0],pRight,SR,cRight,sStar);

    // var 2

    double pLR = 0.5 * (pLeft + pRight + solLeft[0] * cLeft * (sStar - solLeft[1] / solLeft[0]) + \
                                         solRight[0] * cRight * (sStar - solRight[1] / solRight[0]));

    if (sStar >= 0.0)
        return getFStar(solLeft,fluxL,solLeft[0],pLR,SL,cLeft,sStar);

    return getFStar(solRight,fluxR,solLeft[0],pLR,SR,cRight,sStar);


}

