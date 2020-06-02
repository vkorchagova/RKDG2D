#ifndef FLUXVISCOUS_H
#define FLUXVISCOUS_H

#include <algorithm>
#include "Flux.h"

/// 
/// 
/// 

class FluxViscous : public Flux
{

public:

    /// Constructor
    FluxViscous(const Physics& phs) : Flux(phs) {}

    /// Destructor
    ~FluxViscous() {}
    
    /// Evaluate numerical flux through one point
    virtual numvector<double, dimPh> evaluate(
        const numvector<double, dimPh>& solInner, const numvector<double, dimPh>& solOuter,
        const numvector<double, dimGrad>& gradSolInner, const numvector<double, dimGrad>& gradSolOuter
    ) const override;
};


#endif // FLUXVISCOUS_H
