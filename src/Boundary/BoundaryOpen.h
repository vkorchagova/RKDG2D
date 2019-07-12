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
    ~BoundaryOpen() {}
    
    /// Apply boundary condition
    //numvector<double, 5> applyBoundary(const numvector<double, 5>& solLeft = {0.0, 0.0, 0.0, 0.0, 0.0}, const Point& n = Point({0.0,0.0}), int numGP = 0) const override;

    // Get solution in outer side
    virtual numvector<double, dimPh> getSolOuter(
        const numvector<double, dimPh>& solInner, 
        const Point& n = Point({0.0,0.0})
    ) const override;

};

#endif // BOUNDARYOPEN_H
