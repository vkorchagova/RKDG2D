#include "BoundaryNonSlip.h"
#include "Basis.h"


numvector<double, dimPh> BoundaryNonSlip::getSolOuter (
    const numvector<double, dimPh>& solLeft, 
    const Point& n
) const
{    
    return { solLeft[0], -solLeft[1], 0.0, 0.0, solLeft[4] - 0.5 * solLeft[2] * solLeft[2] / solLeft[0] - 0.5 * solLeft[3] * solLeft[3] / solLeft[0] };
    //return { -solLeft[0] + 0.028, -solLeft[1] - 0.28, 0.0, 0.0, solLeft[4] + 0.7};
}

numvector<double, dimGrad> BoundaryNonSlip::getGradSolOuter (
    const numvector<double, dimGrad>& gradSol,
    const Point& n
) const
{
    return gradSol;// {gradSol[0], gradSol[1], gradSol[2], gradSol[3], gradSol[4], gradSol[5], 0, 0};//gradSol;//{0.0, 0.0, 0.0, 0.0, 0.0};
}
