#ifndef LIMITERRIEMANNWENOS_H
#define LIMITERRIEMANNWENOS_H

#include "Limiter.h"

///
/// WENO\_S limiter applied to local characteristic variables component-by-component
///

class LimiterRiemannWENOS : public Limiter
{

public:
    /// Construct by indicator and problem
    LimiterRiemannWENOS(const Indicator& ind, Problem& prb) : Limiter(ind, prb) {}

    /// Limit solution gradients
    ///
    /// @param alpha    vector of solution coeffitients in all cells which should be limited
    virtual void limit(std::vector<numvector<double, 5 * nShapes> >& alpha) override;
};

#endif // LIMITERRIEMANNWENOS_H
