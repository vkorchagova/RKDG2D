#include "BoundaryOpen.h"

numvector<double, dimPh> BoundaryOpen::getSolOuter (
    const numvector<double, dimPh>& solLeft, 
    const Point& n
) const
{
    return solLeft;
}
