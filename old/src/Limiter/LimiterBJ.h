#ifndef LIMITERBJ_H
#define LIMITERBJ_H

#include "Limiter.h"

///
/// WENO\_S limiter applied to conservative variables component-by-component
///

class LimiterBJ : public Limiter
{

public:
    /// Construct by indicator and problem
    LimiterBJ(const Indicator& ind, Problem& prb) : Limiter(ind, prb) {}

    /// Limit solution gradients
    ///
    /// @param alpha    vector of solution coeffitients in all cells which should be limited
    virtual void limit(std::vector<numvector<double, 5 * nShapes> >& alpha) override;
};

#endif // LIMITERBJ_H
