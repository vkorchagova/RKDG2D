#include "BoundaryConstant.h"

numvector<double, 5> BoundaryConstant::applyBoundary(const numvector<double, 5>& solLeft, const Point& n, int numGP ) const
{
    return fixedValue;
}
