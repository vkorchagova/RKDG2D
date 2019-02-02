#include "BoundarySlip.h"
#include "Basis.h"


numvector<double, dimPh> BoundarySlip::getSolOuter (
    const numvector<double, dimPh>& solLeft, 
    const Point& n
) const
{
    numvector<double, dimPh> res = solLeft;
    res[1] *= -1;
    
    return res;
}