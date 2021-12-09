#include "BoundaryOpenTotalPressure.h"

numvector<double, dimPh> BoundaryOpenTotalPressure::getSolOuter (
    const numvector<double, dimPh>& solLeft, 
    const Point& n
) const
{
    // double rhoU2 = solLeft[1]*solLeft[1] + solLeft[2]*solLeft[2];
    // rhoU2 /= solLeft[0];

    // double M = phs.getMachNumber(solLeft);

    // if (M < 1.0)
    // {
    //     double rhoU2 = solLeft[1]*solLeft[1] + solLeft[2]*solLeft[2];
    //     rhoU2 /= solLeft[0];

    //     double energy = phs.e(
    //                     solLeft[0], 
    //                     solLeft[1]/solLeft[0],
    //                     solLeft[2]/solLeft[0], 
    //                     0.0, 
    //                     2.0*pressure - phs.getPressure(solLeft)
    //                 );

    //     if (solLeft[1] < 0)
    //         return { solLeft[0], solLeft[1], solLeft[2], solLeft[3], energy - 0.5*rhoU2 / (phs.cpcv - 1.0) };
    //     else
    //         return { solLeft[0], solLeft[1], solLeft[2], solLeft[3], energy };
    // }

    // return solLeft;

    double rho2U2 = solLeft[1]*solLeft[1] + solLeft[2]*solLeft[2];
    double U2 = rho2U2 / solLeft[0] / solLeft[0];
    double epsIn = solLeft[1] > 0 ? (solLeft[4] - 0.5*rho2U2) / solLeft[0] : phs.R / (phs.cpcv - 1.0) * TFix;

    double M = phs.getMachNumber(solLeft);

    double pressure = M > 1 ? pTotal * pow( 1.0 + 0.5 * (phs.cpcv - 1.0)*M*M, - phs.cpcv / (phs.cpcv - 1.0)) : pTotal;
    double rho = M > 1 ? pressure / (phs.cpcv - 1.0) / epsIn : pressure / ( 0.5*U2 + (phs.cpcv - 1.0) * epsIn);
    pressure = M > 1 ? pressure : pressure - 0.5*U2*rho;

    // double pressure = pTotal * pow( 1.0 + 0.5 * (phs.cpcv - 1.0)*M*M, - phs.cpcv / (phs.cpcv - 1.0));
    // double rho = pressure / (phs.cpcv - 1.0) / epsIn ;    

    double rhoEOut = pressure  / (phs.cpcv - 1.0) + 0.5 * rho * U2;

    return {rho, solLeft[1] / solLeft[0] * rho, solLeft[2] / solLeft[0] * rho, 0.0, rhoEOut};
}
