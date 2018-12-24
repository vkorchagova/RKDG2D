#ifndef LIMITER_H
#define LIMITER_H

//#include "Indicator.h"
#include "numvector.h"
#include "Params.h"
#include <vector>

class Limiter
{

protected:

    //- Discontinuities checker
    //const Indicator& indicator;

    //- Problem
    //Problem& problem;

    //- Number of limitation steps
    static const int nIter = 2;

public:

    //- Construct by indicator and problem
    Limiter() {}
    //Limiter(const Indicator& ind, Problem& prb) : indicator(ind), problem(prb) {}

    //- Limit solution gradients
    void limit(std::vector<numvector<double, dimS>>& alpha) {}

    //virtual void limit(std::vector<numvector<double, 5 * nShapes>>& alpha) = 0;

};

#endif // LIMITER_H
