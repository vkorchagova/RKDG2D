#include "BoundarySine.h"

const double pi = 3.141592653589793;

numvector<double, 5> BoundarySine::applyBoundary(const numvector<double, 5>& sol) const
{
    double sine = a_ * sin(2.0 * pi * f_ * time_.runTime());

    return {sine + u0_[0], sol[1], sol[2], sol[3], sol[4]};


}
