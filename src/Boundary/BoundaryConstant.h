#ifndef BOUNDARYCONSTANT_H
#define BOUNDARYCONSTANT_H

#include "Boundary.h"

/// 
/// Boundary condition: constant value for all conservative variables
/// 

class BoundaryConstant : public Boundary
{
    numvector<double, dimPh> fixedValue;

public:

    /// Constructor
    BoundaryConstant(const Patch& p, const numvector<double, dimPh>& defValue) : Boundary(p), fixedValue(defValue) { type = "constant"; }

    // Get solution in outer side
    virtual numvector<double, dimPh> getSolOuter(
        const numvector<double, dimPh>& solInner, 
        const Point& n = Point({0.0,0.0})
    ) const override;

    virtual numvector<double, dimGrad> getGradSolOuter(
        const numvector<double, dimGrad>& gardSol,
        const Point& n = Point({0.0,0.0})
    ) const override;

};

#endif // BOUNDARYCONSTANT_H
