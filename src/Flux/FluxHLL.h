#ifndef FLUXHLL_H
#define FLUXHLL_H

#include "Flux.h"

/// 
/// HLL numerical flux
/// 

class FluxHLL : public Flux
{
public:

    /// Constructor
    FluxHLL(const Physics& phs) : Flux(phs) {}

    /// Destructor
    ~FluxHLL() {}

    /// Evaluate numerical flux through one point
    virtual numvector<double, dimPh> evaluate(
        const numvector<double, dimPh>& solInner, 
        const numvector<double, dimPh>& solOuter
    ) const override;

};

#endif // FLUXHLL_H
