#ifndef LIMITERFINDIFF_H
#define LIMITERFINDIFF_H

#include "Limiter.h"
#include <vector>

///
/// Simple limiter - just remove solution gradients and make solution piecewise-constant
///

class LimiterFinDiff : public Limiter
{

public:

    /// Construct 
    //LimiterFinDiff(const Indicator& ind, Problem& prb) : Limiter(ind, prb) {}
    LimiterFinDiff(
        const std::vector<std::shared_ptr<Cell>>& cells, 
        const Solution& sln,
        const Physics& phs) : Limiter(cells, sln, phs) {}

    ~LimiterFinDiff() {}

    /// Limit solution gradients
    virtual void limit(std::vector<numvector<double, dimS> >& alpha) override;
};

#endif // LIMITERFINDIFF_H
