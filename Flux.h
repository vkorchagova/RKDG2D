#pragma once

#include "Problem.h"



namespace std
{

class Flux
{

public:
    Problem *problem;

public:
    Flux() {}

    Flux(Problem &prb);
	~Flux();

    virtual vector<double> evaluateHor(numvector<double,2> point, vector<numvector<double,5*nShapes>>& alpha, int iCellUp, int iCellDown) = 0;
    virtual vector<double> evaluateVer(numvector<double,2> point, vector<numvector<double,5*nShapes>>& alpha, int iCellLeft, int iCellRight) = 0;

    //- Get right-hand side
    vector<numvector<double,5*nShapes>> getRHS();

};// end Flux

} // end namespace std

