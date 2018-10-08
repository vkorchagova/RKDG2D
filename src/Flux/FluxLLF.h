#ifndef FLUXLLF_H
#define FLUXLLF_H

#include <algorithm>
#include "Flux.h"

///
/// Local Lax --- Fridriechs Riemann solver
///

class FluxLLF : public Flux
{

public:

    /// Constructor
    FluxLLF(const Problem& prb) : Flux(prb) {}

    /// Destructor
    ~FluxLLF() {}

    /// Evaluate numerical flux through one point
    ///
    /// @param solInner    solution at one side of edge
    /// @param solOuter    solution at other side of edge
    /// @param n           normal direction to edge
    ///
    /// @return    values of numerical fluxes of all conservative variables
    virtual numvector<double, 5> evaluate(const numvector<double, 5>& solInner, const numvector<double, 5>& solOuter, const Point& n) const override;
};


#endif // FLUXLLF_H
