#ifndef FLUXHLLC_H
#define FLUXHLLC_H

#include "Flux.h"

///
/// HLLC Riemann solver
///

class FluxHLLC : public Flux
{
    /// Calculate Star Region velocity
    numvector<double, 5> getUStar (const numvector<double, 5>& sol, double pK, double SK, double cK, double SStar) const;

    /// Calculate auxiliary variable used in Star Region velocity computation
    numvector<double, 5> getFStar(const numvector<double, 5>& sol, const numvector<double, 5>& fK, double rhoL, double pK, double SK, double cK, double SStar) const;

public:

    /// Constructor
    FluxHLLC(const Problem& prb) : Flux(prb) {}

    /// Destructor
    ~FluxHLLC() {}

    /// Evaluate numerical flux through one point
    ///
    /// @param solInner    solution at one side of edge
    /// @param solOuter    solution at other side of edge
    /// @param n           normal direction to edge
    ///
    /// @return    values of numerical fluxes of all conservative variables
    virtual numvector<double, 5> evaluate(const numvector<double, 5>& solLeft, const numvector<double, 5>& solRight, const Point& n) const override;

};

#endif // FLUXHLLC_H
