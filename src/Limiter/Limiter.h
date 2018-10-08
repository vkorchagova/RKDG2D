#ifndef LIMITER_H
#define LIMITER_H

#include "Indicator.h"
#include "Problem.h"
#include <vector>

///
/// Abstract class for troubled cell indicators
///

class Limiter
{

protected:

    /// Constant reference to indicator of discontinuities
    const Indicator& indicator;

    /// Reference to problem
    Problem& problem;

    /// Number of limitation steps
    static const int nIter = 2;

public:

    /// Construct by indicator and problem
    Limiter(const Indicator& ind, Problem& prb) : indicator(ind), problem(prb) {}

    /// Limit solution gradients
    ///
    /// @param alpha    vector of solution coeffitients in all cells which should be limited
    virtual void limit(std::vector<numvector<double, 5 * nShapes>>& alpha) = 0;

};

#endif // LIMITER_H
