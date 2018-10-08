#ifndef LIMITERWENOS_H
#define LIMITERWENOS_H

#include "Limiter.h"

///
/// WENO\_S limiter applied to conservative variables component-by-component
///

class LimiterWENOS : public Limiter
{

public:
    /// Construct by indicator and problem
    LimiterWENOS(const Indicator& ind, Problem& prb) : Limiter(ind, prb) {}

    /// Limit solution gradients
    ///
    /// @param alpha    vector of solution coeffitients in all cells which should be limited
    virtual void limit(std::vector<numvector<double, 5 * nShapes> >& alpha) override;
};

#endif // LIMITERWENOS_H
