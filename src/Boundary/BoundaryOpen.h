#ifndef BOUNDARYOPEN_H
#define BOUNDARYOPEN_H

#include "Boundary.h"

/// Open boundary condition
///
/// Provides the similar values of conservative variables inside and outside the boundary cell
///
/// Conservative variables are:
/// -   density;
/// -   normal component of velocity;
/// -   tangential component of velocity;
/// -   z-component of velocity (0);
/// -   volumetric total energy.

class BoundaryOpen : public Boundary
{
public:

    /// Default constructor
    BoundaryOpen() : Boundary() { type = "open";}

    /// Apply boundary condition
    ///
    /// @param  sol solution near the point n inside boundary cell
    /// @param  n   Normal direction related to boundary edge (not used)
    ///
    /// @return the value of solution equals to sol
    numvector<double, 5> applyBoundary(const numvector<double, 5>& solLeft = {0.0, 0.0, 0.0, 0.0, 0.0}, const Point& n = Point({0.0,0.0})) const override;

};

#endif // BOUNDARYOPEN_H
