#ifndef FLUX_H
#define FLUX_H

#include "Problem.h"
#include "defs.h"

///
/// Abstract class for Riemann solvers
///

class Flux
{

public:

    /// Reference to problem
    const Problem& problem;

public:

    /// Constructor
    Flux(const Problem& prb) : problem(prb) {};

    /// Destructor
    virtual ~Flux() {}

    /// Evaluate numerical flux through one point
    ///
    /// @param solInner    solution at one side of edge
    /// @param solOuter    solution at other side of edge
    /// @param n           normal direction to edge
    ///
    /// @return    values of numerical fluxes of all conservative variables
    virtual numvector<double, 5> evaluate(const numvector<double, 5>& solInner, const numvector<double, 5>& solOuter, const Point& n) const = 0;

};// end Flux

#endif // FLUX_H

