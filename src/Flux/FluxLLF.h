/// -----------------------------------
/// Local Lax-Fridriechs numerical flux
/// -----------------------------------



#ifndef FLUXLLF_H
#define FLUXLLF_H

#include "Flux.h"

class FluxLLF : public Flux
{

public:

    //- Default constructor
    FluxLLF() : Flux() {}

    //- Construct with problem
    FluxLLF(Problem &prb) : Flux(prb) {}

    //- Destructor
    ~FluxLLF() {}
    
    //- Overloaded "=" operator
    FluxLLF& operator= (const Flux& flx)  { problem = flx.problem; return *this; }

    //- Evaluate numerical flux through one point
    virtual numvector<double, 5> evaluate(const numvector<double, 5>& solInner, const numvector<double, 5>& solOuter, const Point& n) const override;
};


#endif // FLUXLLF_H
