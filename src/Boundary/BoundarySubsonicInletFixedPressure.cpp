#include "BoundarySubsonicInletFixedPressure.h"
#include "compService.h"

numvector<double, dimPh> BoundarySubsonicInletFixedPressure::getSolOuter (
    const numvector<double, dimPh>& solLeft, 
    const Point& n
) const
{
    double rhoU2 = solLeft[1]*solLeft[1] + solLeft[2]*solLeft[2];
    rhoU2 /= solLeft[0];

    double MTTotSq = rhoU2 / solLeft[0] / phs.cpcv / phs.R / TTot;

    double MSq = MTTotSq / (1.0 - 0.5 * (phs.cpcv - 1.0) * MTTotSq);

    // fixedValue
    double p = pFix;

    double T = TTot / ( 1.0 + 0.5 * (phs.cpcv - 1.0)*MSq );

    double rho = p / phs.R / T;
    
    return { rho, 
             solLeft[1] / solLeft[0] * rho, 
             solLeft[2] / solLeft[0] * rho, 
             0.0, 
             p / (phs.cpcv - 1.0) + 0.5 * rhoU2 / solLeft[0] * rho};
}
