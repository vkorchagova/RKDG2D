#ifndef LIMITERFINDIFF_H
#define LIMITERFINDIFF_H

#include "Limiter.h"
#include <vector>

///
/// FinDiff limiter
///
/// Delete gradients of solution in troubled cells
///


class LimiterFinDiff : public Limiter
{

public:

    /// Construct by indicator and problem
    LimiterFinDiff(const Indicator& ind, Problem& prb) : Limiter(ind, prb) {}

    /// Limit solution gradients
    ///
    /// @param alpha    vector of solution coeffitients in all cells which should be limited
    virtual void limit(std::vector<numvector<double, 5 * nShapes> >& alpha) override;
};

#endif // LIMITERFINDIFF_H
