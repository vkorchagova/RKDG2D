#include "BoundaryOpen.h"

numvector<double, 5> BoundaryOpen::applyBoundary(const numvector<double, 5>& solLeft, const Point& n) const
{
    return { solLeft[0], solLeft[1], solLeft[2], solLeft[3], solLeft[4]};
}
