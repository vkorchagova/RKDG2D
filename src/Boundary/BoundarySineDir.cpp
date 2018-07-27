#include "BoundarySineDir.h"
#include "defs.h"

const double pi = 3.141592653589793;

numvector<double, 5> BoundarySineDir::applyBoundary(const numvector<double, 5>& sol, const Point& n) const
{
    double sine = a_ * sin(2.0 * pi * f_ * time_.runTime());

    return inverseRotate({sol[0], sine + u0_[1], u0_[2], sol[3], sol[4]}, n);
}
