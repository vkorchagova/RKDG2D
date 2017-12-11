#ifndef FLUXLLF_H
#define FLUXLLF_H

#include "Flux.h"

class FluxLLF : public Flux
{
public:
    FluxLLF() : Flux() {}
    FluxLLF(Problem &prb) : Flux(prb) {}
    ~FluxLLF() {}

    virtual numvector<double,5> evaluate(numvector<double,5>& solLeft, numvector<double,5>& solRight);

    //numvector<double,5> evaluateHor(numvector<double,2> point, const vector<numvector<double,5*nShapes>>& alpha, int iCellUp, int iCellDown);
    //numvector<double,5> evaluateVer(numvector<double,2> point, const vector<numvector<double,5*nShapes>>& alpha, int iCellLeft, int iCellRight);

    //- Get right-hand side
    //std::vector<numvector<double,5*nShapes>> getRHS(const std::vector<numvector<double,5*nShapes>>& alpha);
};


#endif // FLUXLLF_H
