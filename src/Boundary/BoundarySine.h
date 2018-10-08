#ifndef BOUNDARYSINE_H
#define BOUNDARYSINE_H

#include "numvector.h"
#include "TimeClass.h"
#include "Boundary.h"

class Problem;

/// Fixed value boundary condition (pulsating velocity in the normal direction to boundary edge)
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
/// -   normal component of velocity;
/// -   tangential component of velocity;
/// -   z-component of velocity (0);
/// -   volumetric total energy.

class BoundarySine : public Boundary
{

protected:

    /// Amplitude
    double a_;

    /// Frequency
    double f_;

    /// Reference value of all conservative variables
    numvector<double, 5> u0_;

    /// Reference to time object
    const Time& time_;

    /// Reference to problem
    const Problem& problem_;

public:

    /// Constructor
    BoundarySine(double a, double f, const Time& t, const Problem& prb, const numvector<double,5>& u0 = {1,0,0,0,1}) : Boundary(), a_(a), f_(f), u0_(u0), time_(t), problem_(prb) { type = "sine";}

    // Set reference value
    //void setRefValue(const numvector<double,5>& u0) {u0_ = u0;}

    /// Apply boundary condition
    ///
    /// @param  sol solution near the point n inside boundary cell 
    /// @param  n   Normal direction related to boundary edge (not used)
    ///
    /// @return the calculated value of solution 
    virtual numvector<double, 5> applyBoundary(const numvector<double, 5>& solLeft = {0.0, 0.0, 0.0, 0.0, 0.0}, const Point& n = Point({0.0,0.0})) const override;

};

#endif // BOUNDARYSINE_H
