#ifndef FLUXLLF_H
#define FLUXLLF_H

#include "Flux.h"

class FluxLLF : public Flux
{
public:
    FluxLLF() : Flux() {}
    FluxLLF(Problem &prb) : Flux(prb) {}
    ~FluxLLF() {}
    
    FluxLLF& operator= (const Flux& flx)  { problem = flx.problem; return *this; }

    virtual numvector<double,5> evaluate(const numvector<double,5>& solLeft, const numvector<double,5>& solRight) override;

    //numvector<double,5> evaluateHor(numvector<double,2> point, const vector<numvector<double,5*nShapes>>& alpha, int iCellUp, int iCellDown);
    //numvector<double,5> evaluateVer(numvector<double,2> point, const vector<numvector<double,5*nShapes>>& alpha, int iCellLeft, int iCellRight);

    //- Get right-hand side
    //std::vector<numvector<double,5*nShapes>> getRHS(const std::vector<numvector<double,5*nShapes>>& alpha);
};


#endif // FLUXLLF_H
