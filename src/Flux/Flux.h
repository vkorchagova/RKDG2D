#ifndef FLUX_H
#define FLUX_H

#include "Problem.h"

class Flux
{

public:
    Problem *problem;

public:
    Flux() {}

    Flux(Problem &prb);
    ~Flux() {}

    virtual numvector<double,5> evaluate(numvector<double,5>& , numvector<double,5>& ) {}

//    virtual numvector<double,5> evaluateHor(numvector<double,2> point, const vector<numvector<double,5*nShapes>>& alpha, int iCellUp, int iCellDown) = 0;
//    virtual numvector<double,5> evaluateVer(numvector<double,2> point, const vector<numvector<double,5*nShapes>>& alpha, int iCellLeft, int iCellRight) = 0;

    //- Get right-hand side
    //virtual vector<numvector<double,5*nShapes>> getRHS(const vector<numvector<double,5*nShapes>>& alpha) = 0;

};// end Flux

#endif // FLUX_H

