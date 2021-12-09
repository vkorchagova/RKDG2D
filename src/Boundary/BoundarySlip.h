#ifndef BOUNDARYSLIP_H
#define BOUNDARYSLIP_H

#include "Boundary.h"

/// 
/// Slip boundary condition
/// 

class BoundarySlip : public Boundary
{

public:

    /// Constructor
    BoundarySlip(const Patch& p, const Physics& _phs) : Boundary(p, _phs) { type = "slip"; }

    // Get solution in outer side
    virtual numvector<double, dimPh> getSolOuter(
        const numvector<double, dimPh>& solInner, 
        const Point& n = Point({0.0,0.0})
    ) const override;
};

#endif // BOUNDARYSLIP_H

