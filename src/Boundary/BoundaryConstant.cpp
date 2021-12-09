#include "BoundaryConstant.h"
#include "compService.h"

numvector<double, dimPh> BoundaryConstant::getSolOuter (
    const numvector<double, dimPh>& solLeft, 
    const Point& n
) const
{
    return 2.0*fixedValue - solLeft;
}
