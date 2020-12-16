#ifndef BOUNDARYINLET_H
#define BOUNDARYINLET_H

#include "Boundary.h"

///
/// Inlet boundary condition
///

class BoundaryInlet : public Boundary
{

    numvector<double, dimPh> fixedValue;

public:

    /// Constructor
    BoundaryInlet(const Patch& p, const numvector<double, dimPh>& defValue) : Boundary(p), fixedValue(defValue) { type = "inlet"; }

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

#endif // BOUNDARYINLET_H

