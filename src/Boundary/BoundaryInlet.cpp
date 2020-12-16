#include "BoundaryInlet.h"
#include "Basis.h"


numvector<double, dimPh> BoundaryInlet::getSolOuter (
    const numvector<double, dimPh>& solLeft,
    const Point& n
) const
{
    //return { -solLeft[0] + 0.028, -solLeft[1] - 0.28, 0.0, 0.0, solLeft[4] + 0.7};
    return { -solLeft[0] + 2.0 * fixedValue[0], -solLeft[1] + 2.0 * fixedValue[1], -solLeft[2] + 2.0 * fixedValue[2], \
                -solLeft[3] + 2.0 * fixedValue[3],\
                -solLeft[4] + 2.0 * fixedValue[4] };
                //(solLeft[4] - 0.5 * (solLeft[1] * solLeft[1] + solLeft[2] * solLeft[2]) / solLeft[0]) \
                + 2.0 * 0.5 * (fixedValue[1] * fixedValue[1] + fixedValue[2] * fixedValue[2]) / fixedValue[0]};
}

numvector<double, dimGrad> BoundaryInlet::getGradSolOuter (
    const numvector<double, dimGrad>& gradSol,
    const Point& n
) const
{
    return  gradSol;//{0.0, 0.0, 0.0, 0.0, 0.0};
}


