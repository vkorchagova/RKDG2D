#ifndef FLUXLLF_H
#define FLUXLLF_H

#include <algorithm>
#include "Flux.h"

/// 
/// Local Lax-Fridriechs numerical flux
/// 

class FluxLLF : public Flux
{

public:

    /// Constructor
    FluxLLF(const Physics& phs) : Flux(phs) {}

    /// Destructor
    ~FluxLLF() {}
    
    /// Evaluate numerical flux through one point
    virtual numvector<double, dimPh> evaluate(const numvector<double, dimPh>& solInner, 
											  const numvector<double, dimPh>& solOuter) const override;
};


#endif // FLUXLLF_H
