#include "BoundaryConstant.h"
#include "compService.h"

numvector<double, dimPh> BoundaryConstant::getSolOuter (
    const numvector<double, dimPh>& solLeft, 
    const Point& n
) const
{
    return fixedValue;
}

numvector<double, dimGrad> BoundaryConstant::getGradSolOuter(
    const numvector<double, dimGrad>& gradSol,
    const Point& n
) const
{
    return 0.0; // !!!! grad fixedValue
}
