#ifndef BOUNDARYOPEN_H
#define BOUNDARYOPEN_H

#include "Boundary.h"

/// 
/// Open boundary condition
/// 

class BoundaryOpen : public Boundary
{
public:

    /// Default constructor
    BoundaryOpen(const Patch& p) : Boundary(p) { type = "open";}

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

#endif // BOUNDARYOPEN_H
