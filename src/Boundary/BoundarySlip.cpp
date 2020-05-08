#include "BoundarySlip.h"
#include "Basis.h"


numvector<double, dimPh> BoundarySlip::getSolOuter (
    const numvector<double, dimPh>& solLeft, 
    const Point& n
) const
{    
    return { solLeft[0], -solLeft[1], solLeft[2], solLeft[3], solLeft[4]};
}

numvector<double, dimGrad> BoundarySlip::getGradSolOuter (
    const numvector<double, dimGrad>& gradSol,
    const Point& n
) const
{
    return  gradSol;
}
