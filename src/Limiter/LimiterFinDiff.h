#ifndef LIMITERFINDIFF_H
#define LIMITERFINDIFF_H

#include "Limiter.h"
#include <vector>

///
/// Simple limiter - just remove solution gradients and make solution piecewise-constant
///

class LimiterFinDiff : public Limiter
{

protected:
    
    /// Limit solution gradients
    virtual numvector<double, dimS> limitation(const std::vector<std::shared_ptr<Cell>>& stencil) override;

    /// Choose stencil for defined cell
    virtual std::vector<std::shared_ptr<Cell>> getStencilFor(const std::shared_ptr<Cell>& cell) override;

public:

    /// Construct 
    //LimiterFinDiff(const Indicator& ind, Problem& prb) : Limiter(ind, prb) {}
    LimiterFinDiff(
        const Mesh& msh,
        Solution& sln,
        const Physics& phs,
        const Indicator& ind,
        Buffers& buf) : Limiter(msh, sln, phs, ind, buf) {};

    ~LimiterFinDiff() {}
};

#endif // LIMITERFINDIFF_H
