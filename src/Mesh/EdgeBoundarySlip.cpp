#include "EdgeBoundarySlip.h"

numvector<double, 5> EdgeBoundarySlip::applyBoundary(const numvector<double, 5>& solLeft) const
{
    return { solLeft[0], solLeft[1], solLeft[2], solLeft[3], solLeft[4]};
}
