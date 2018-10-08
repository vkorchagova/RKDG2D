#ifndef LIMITERMUSCL_H
#define LIMITERMUSCL_H

#include "Limiter.h"

///
/// MUSCL limiter
///

class LimiterMUSCL : public Limiter
{
public:
    /// Construct by indicator and problem
    LimiterMUSCL(const Indicator& ind, Problem& prb) : Limiter(ind, prb) {}

    /// Limit solution gradients
    ///
    /// @param alpha    vector of solution coeffitients in all cells which should be limited
    virtual void limit(std::vector<numvector<double, 5 * nShapes> >& alpha) override;
};

#endif // LIMITERMINMOD_H
