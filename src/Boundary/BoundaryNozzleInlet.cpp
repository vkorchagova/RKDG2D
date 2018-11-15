#include "BoundaryNozzleInlet.h"

numvector<double, 5> BoundaryNozzleInlet::applyBoundary(const numvector<double, 5>& solLeft, const Point& n) const
{
    double rhoSol2 = solLeft[0] * solLeft[0];
    double magU2 = (solLeft[1] * solLeft[1] + solLeft[2] * solLeft[2]) / rhoSol2;
    
    double R = R0 / M;
    
    
    double rho = pTot / (R * T + magU2 * 0.25 * (2.0 * problem.cpcv - 3) / (problem.cpcv - 1.0));
    
    //double p = pTot - 0.5 * rho * magU2;
    
    double e = rho * R * T / (problem.cpcv - 1.0);
    
    double corrMom = rho / solLeft[0];
    
    
    return { rho, -corrMom * solLeft[1], corrMom * solLeft[2], 0.0, e};
}
