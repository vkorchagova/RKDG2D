#ifndef FLUXLLF_H
#define FLUXLLF_H

#include "Flux.h"

namespace std
{

class FluxLLF : public Flux
{
public:
    FluxLLF();
    ~FluxLLF();

    vector<double> evaluateHor(numvector<double,2> point, vector<numvector<double,5*nShapes>>& alpha, int iCellUp, int iCellDown);
    vector<double> evaluateVer(numvector<double,2> point, vector<numvector<double,5*nShapes>>& alpha, int iCellLeft, int iCellRight);
};

}

#endif // FLUXLLF_H
