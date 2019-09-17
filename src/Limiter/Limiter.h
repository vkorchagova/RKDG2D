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

///
/// Abstract class for limiter of solution
///

class Limiter
{

protected:

    /// Discontinuities checker
    //const Indicator& indicator;

    /// Constant reference to mesh
    const std::vector<std::shared_ptr<Cell>>& cells;

    /// Reference to solution
    const Solution& solution;

    /// Problem
    const Physics& physics;

    /// Number of limitation steps
    static const int nIter = 2;

    /// List of numbers of troubled cells
    std::vector<int> troubledCells;

    /// Limited solution
    std::vector<numvector<double, dimS>> alphaNew;

    /// Last hope limiter
    void lastHope(std::vector<numvector<double, dimS> >& alpha);

public:

    /// Construct
    Limiter(const std::vector<std::shared_ptr<Cell>>& cells, 
        const Solution& sln,
        const Physics& phs);

    /// Destructor
    virtual ~Limiter() {}

    /// Limit solution gradients
    virtual void limit(std::vector<numvector<double, dimS> >& alpha) = 0;

    

};

#endif // LIMITER_H
