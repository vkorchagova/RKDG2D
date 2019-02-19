#ifndef BOUNDARYSLIP_H
#define BOUNDARYSLIP_H

#include "Boundary.h"

class BoundarySlip : public Boundary
{

public:

    /// Constructor
    BoundarySlip(const Patch& p) : Boundary(p) { type = "slip"; }
    ~BoundarySlip() {}

    /// Apply boundary condition
    //virtual void applyBoundary(std::vector<numvector<double, dimS>>& coeffs) const override;

    virtual numvector<double, dimPh> getSolOuter(
        const numvector<double, dimPh>& solInner, 
        const Point& n = Point({0.0,0.0})
    ) const override;
};

#endif // BOUNDARYSLIP_H

