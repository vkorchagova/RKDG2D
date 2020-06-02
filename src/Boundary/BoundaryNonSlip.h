#ifndef BOUNDARYNONSLIP_H
#define BOUNDARYNONSLIP_H

#include "Boundary.h"

/// 
/// Non-Slip boundary condition
/// 

class BoundaryNonSlip : public Boundary
{

public:

    /// Constructor
    BoundaryNonSlip(const Patch& p) : Boundary(p) { type = "non-slip"; }

    // Get solution in outer side
    virtual numvector<double, dimPh> getSolOuter(
        const numvector<double, dimPh>& solInner, 
        const Point& n = Point({0.0,0.0})
    ) const override;

    virtual numvector<double, dimGrad> getGradSolOuter(
        const numvector<double, dimGrad>& GradSol,
        const Point& n = Point({0.0,0.0})
    ) const override;
};

#endif // BOUNDARYNONSLIP_H

