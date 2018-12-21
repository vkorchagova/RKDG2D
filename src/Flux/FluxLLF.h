/// -----------------------------------
/// Local Lax-Fridriechs numerical flux
/// -----------------------------------


#ifndef FLUXLLF_H
#define FLUXLLF_H


#include <algorithm>
#include "Flux.h"

class FluxLLF : public Flux
{

public:

    //- Default constructor
    // FluxLLF() : Flux() {}

    //- Construct with problem
    FluxLLF(const Physics& phs) : Flux(phs) {}

    //- Destructor
    ~FluxLLF() {}
    
    //- Evaluate numerical flux through one point
    virtual numvector<double, dimPh> evaluate(const numvector<double, dimPh>& solInner, 
											  const numvector<double, dimPh>& solOuter, const Point& n) const override;
};


#endif // FLUXLLF_H
