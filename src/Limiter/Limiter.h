#ifndef LIMITER_H
#define LIMITER_H

//#include "Indicator.h"
#include "numvector.h"
#include "Params.h"
#include "Cell.h"
#include "Solution.h"
#include "Physics.h"
#include <vector>
#include <memory>

class Limiter
{

protected:

    //- Discontinuities checker
    //const Indicator& indicator;

    //- Constant reference to mesh
    const std::vector<std::shared_ptr<Cell>>& cells;

    //- Reference to solution
    const Solution& solution;

    //- Problem
    const Physics& physics;

    //- Number of limitation steps
    static const int nIter = 2;

public:

    //- Construct
    Limiter(
        const std::vector<std::shared_ptr<Cell>>& cells, 
        const Solution& sln,
        const Physics& phs) : cells(cells), solution(sln), physics(phs) {}

    //- Destructor
    virtual ~Limiter() {}
    //Limiter(const Indicator& ind, Problem& prb) : indicator(ind), problem(prb) {}

    //- Limit solution gradients
    //void limit(std::vector<numvector<double, dimS>>& alpha) {}

    virtual void limit(std::vector<numvector<double, dimS> >& alpha) = 0;

    //- Last hope limiter
    void lastHope(std::vector<numvector<double, dimS> >& alpha);

};

#endif // LIMITER_H
