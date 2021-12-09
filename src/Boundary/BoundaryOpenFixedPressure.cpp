#include "BoundaryOpenFixedPressure.h"

numvector<double, dimPh> BoundaryOpenFixedPressure::getSolOuter (
    const numvector<double, dimPh>& solLeft, 
    const Point& n
) const
{
    double rho2U2 = solLeft[1]*solLeft[1] + solLeft[2]*solLeft[2];
    double U2 = rho2U2 / solLeft[0] / solLeft[0];
    double epsIn = (solLeft[4] - 0.5*rho2U2) / solLeft[0];

    double rho = pressure / (phs.cpcv - 1.0) / epsIn;
    double rhoEOut = pressure  / (phs.cpcv - 1.0) + 0.5 * rho * U2;

    return {rho, solLeft[1] / solLeft[0] * rho, solLeft[2] / solLeft[0] * rho, 0.0, rhoEOut};
}
