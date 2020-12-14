#ifndef FLUXPH_H
#define FLUXPH_H

#include <algorithm>
#include "Flux.h"

/// 
/// Physics numerical flux
/// 

class FluxPh : public Flux
{

public:

    /// Constructor
    FluxPh(const Physics& phs) : Flux(phs) {}

    /// Destructor
    ~FluxPh() {}
    
    /// Evaluate numerical flux through one point
    virtual numvector<double, dimPh> evaluate(
        const numvector<double, dimPh>& solInner, const numvector<double, dimPh>& solOuter,
        const numvector<double, dimGrad>& gradSolInner, const numvector<double, dimGrad>& gradSolOuter
    ) const override;
};


#endif // FLUXPH_H
