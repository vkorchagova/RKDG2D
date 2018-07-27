#ifndef BOUNDARYSLIP_H
#define BOUNDARYSLIP_H

#include "Boundary.h"

class BoundarySlip : public Boundary
{
public:
    BoundarySlip() : Boundary() { type = "slip";}

    //- Apply boundary condition
    numvector<double, 5> applyBoundary(const numvector<double, 5>& solLeft = {0.0, 0.0, 0.0, 0.0, 0.0}, const Point& n = Point({0.0,0.0})) const override;
};

#endif // BOUNDARYSLIP_H
