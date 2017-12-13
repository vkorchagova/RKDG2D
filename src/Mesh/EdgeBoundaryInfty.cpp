#include "EdgeBoundaryInfty.h"
#include "Cell.h"

numvector<double, 5> EdgeBoundaryInfty::applyBoundary(const numvector<double, 5>& solLeft) const
{
    return infty;
}
