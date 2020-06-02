#include "BoundaryOpen.h"

numvector<double, dimPh> BoundaryOpen::getSolOuter (
    const numvector<double, dimPh>& solLeft, 
    const Point& n
) const
{
    return solLeft;
}

numvector<double, dimGrad> BoundaryOpen::getGradSolOuter (
    const numvector<double, dimGrad>& gradSolLeft,
    const Point& n
) const
{
    return gradSolLeft;
}
