#ifndef FLUXLLF_H
#define FLUXLLF_H

#include "Flux.h"

namespace std
{

class FluxLLF : public Flux
{
public:
    FluxLLF();
    FluxLLF(Problem &prb);
    ~FluxLLF();

    numvector<double,5> evaluateHor(numvector<double,2> point, const vector<numvector<double,5*nShapes>>& alpha, int iCellUp, int iCellDown);
    numvector<double,5> evaluateVer(numvector<double,2> point, const vector<numvector<double,5*nShapes>>& alpha, int iCellLeft, int iCellRight);

    //- Get right-hand side
    vector<numvector<double,5*nShapes>> getRHS(const vector<numvector<double,5*nShapes>>& alpha);
};

}

#endif // FLUXLLF_H
