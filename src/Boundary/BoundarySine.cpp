#include "BoundarySine.h"

numvector<double, 5> BoundarySine::applyBoundary(const numvector<double, 5>& solLeft) const
{
    double sine = a_ * sin(f_ * time_.runTime());

    return 0;
}
