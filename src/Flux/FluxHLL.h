#ifndef FLUXHLL_H
#define FLUXHLL_H

#include "Flux.h"

class FluxHLL : public Flux
{
public:
    //- Default constructor
    // FluxHLL() : Flux() {}

    //- Construct with problem
    FluxHLL(const Problem& prb) : Flux(prb) {}

    //- Destructor
    ~FluxHLL() {}

    //- Evaluate numerical flux through one point
    virtual numvector<double, 5> evaluate(const numvector<double, 5>& solInner, const numvector<double, 5>& solOuter, const Point& n) const override;

};

#endif // FLUXHLL_H
