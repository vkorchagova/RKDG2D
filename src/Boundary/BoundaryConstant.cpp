#include "BoundaryConstant.h"

numvector<double, 5> BoundaryConstant::applyBoundary(const numvector<double, 5>& solLeft) const
{
    return fixedValue;
}
