#include "BoundaryNonSlip.h"
#include "Basis.h"


numvector<double, dimPh> BoundaryNonSlip::getSolOuter (
    const numvector<double, dimPh>& solLeft, 
    const Point& n
) const
{    
    return { solLeft[0], -solLeft[1], 0.0, 0.0, solLeft[4] - 0.5 * solLeft[2] * solLeft[2] / solLeft[0] - 0.5 * solLeft[3] * solLeft[3] / solLeft[0] };
}

numvector<double, dimGrad> BoundaryNonSlip::getGradSolOuter (
    const numvector<double, dimGrad>& gradSol,
    const Point& n
) const
{
    return  gradSol;//{0.0, 0.0, 0.0, 0.0, 0.0};
}
