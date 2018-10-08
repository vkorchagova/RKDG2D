#ifndef BOUNDARYSLIP_H
#define BOUNDARYSLIP_H

#include "Boundary.h"

/// Zero gradient boundary condition
///
/// Provides zero flow in the normal direction to boundary edge
/// 
/// Conservative variables are:
/// -   density;
/// -   normal component of velocity;
/// -   tangential component of velocity;
/// -   z-component of velocity (0);
/// -   volumetric total energy.

class BoundarySlip : public Boundary
{

public:

    /// Constructor
    BoundarySlip() : Boundary() { type = "slip";}

    /// Apply boundary condition
    ///
    /// @param  sol solution near the point n inside boundary cell
    /// @param  n   Normal direction related to boundary edge (not used)
    ///
    /// @return the calculated value of solution 
    numvector<double, 5> applyBoundary(const numvector<double, 5>& solLeft = {0.0, 0.0, 0.0, 0.0, 0.0}, const Point& n = Point({0.0,0.0})) const override;
};

#endif // BOUNDARYSLIP_H
