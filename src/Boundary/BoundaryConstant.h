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
    ~BoundaryConstant() {};

    /// Apply boundary condition
   // numvector<double, dimPh> applyBoundary(const numvector<double, dimPh>& solLeft = {0.0, 0.0, 0.0, 0.0, 0.0}, const Point& n = Point({0.0,0.0}), int numGP = 0) const override;

    // Get solution in outer side
    virtual numvector<double, dimPh> getSolOuter(
        const numvector<double, dimPh>& solInner, 
        const Point& n = Point({0.0,0.0})
    ) const override;

};

#endif // BOUNDARYCONSTANT_H
