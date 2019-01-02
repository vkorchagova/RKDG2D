#ifndef BOUNDARYSLIP_H
#define BOUNDARYSLIP_H

#include "Boundary.h"

class BoundarySlip : public Boundary
{

public:

    //- Constructor
    BoundarySlip(const Patch& p) : Boundary(p) { type = "slip"; }

    //- Apply boundary condition
    virtual void applyBoundary(std::vector<numvector<double, dimS>>& coeffs) const override;
};

#endif // BOUNDARYSLIP_H

