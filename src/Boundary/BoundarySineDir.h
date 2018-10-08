#ifndef BOUNDARYSINEDIR_H
#define BOUNDARYSINEDIR_H

#include "BoundarySine.h"

/// Fixed value boundary condition (pulsating velocity in the defined direction)
///
/// Provides the pulsating value of velocity on the boundary:
///     \f[
///         u = u_0 + a \sin(2 \pi f t),
///     \f]
/// where a is amplitude, f is frequency, t is time, \f$ u_0 \f$ is reference value of velocity.
///
/// Volumetric total energy held constant.
/// 
/// Conservative variables are:
/// -   density;
/// -   x-component of velocity;
/// -   y-component of velocity;
/// -   z-component of velocity (0);
/// -   volumetric total energy.

class BoundarySineDir : public BoundarySine
{
public:
    
    /// Constructor
    BoundarySineDir(double a, double f, const Time& t, const Problem& prb, const numvector<double,5>& u0 = {1,0,0,0,1}) : BoundarySine(a, f, t, prb, u0) { type = "sineDir";}

    /// Apply boundary condition
    ///
    /// @param  sol solution near the point n inside boundary cell 
    /// @param  n   Normal direction related to boundary edge
    ///
    /// @return the calculated value of solution
    virtual numvector<double, 5> applyBoundary(const numvector<double, 5>& solLeft = {0.0, 0.0, 0.0, 0.0, 0.0}, const Point& n = Point({1.0,0.0})) const override;

};

#endif // BOUNDARYSINEDIR_H
