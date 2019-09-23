#ifndef LIMITER_H
#define LIMITER_H

//#include "Indicator.h"
#include "numvector.h"
#include "Params.h"
#include "Cell.h"
#include "Solution.h"
#include "Physics.h"
#include "Indicator.h"
#include <vector>
#include <memory>

///
/// Abstract class for limiter of solution
///

class Limiter
{
private:

    /// List of numbers of troubled cells
    std::vector<int> troubledCells;

    /// Limited solution
    std::vector<numvector<double, dimS>> newSOL;

    /// Last hope limiter
    void lastHope(std::vector<numvector<double, dimS> >& alpha);

protected:

    /// Discontinuities checker
    const Indicator& indicator;

    /// Problem
    const Physics& physics;

    /// Constant reference to mesh
    const Mesh& mesh;

    /// Reference to solution
    Solution& solution;

    /// Max possible stencil size
    static const int maxPossibleStencilSize = 10;

    /// Limitation algorithm for solution in defined cell
    virtual numvector<double, dimS> limitation(const std::vector<std::shared_ptr<Cell>>& stencil ) = 0;

    /// Update stencil
    virtual std::vector<std::shared_ptr<Cell>> getStencilFor(const std::shared_ptr<Cell>& cell) = 0;

public:

    /// Construct
    Limiter(const Mesh& mesh, 
        Solution& sln,
        const Physics& phs,
        const Indicator& ind);

    /// Destructor
    virtual ~Limiter() {}

    /// Limitation procedure for solution in the entire flow domain
    //void limit(std::vector<numvector<double, dimS> >& alpha);
    void limitSolution();
};

#endif // LIMITER_H
