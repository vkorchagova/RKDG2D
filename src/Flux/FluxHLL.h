#ifndef FLUXHLL_H
#define FLUXHLL_H

#include "Flux.h"

///
/// HLL Riemann solver
///

class FluxHLL : public Flux
{
public:

    /// Construct with problem
    FluxHLL(const Problem& prb) : Flux(prb) {}

    /// Destructor
    ~FluxHLL() {}

    /// Evaluate numerical flux through one point
    ///
    /// @param solInner    solution at one side of edge
    /// @param solOuter    solution at other side of edge
    /// @param n           normal direction to edge
    ///
    /// @return    values of numerical fluxes of all conservative variables
    virtual numvector<double, 5> evaluate(const numvector<double, 5>& solInner, const numvector<double, 5>& solOuter, const Point& n) const override;

};

#endif // FLUXHLL_H
