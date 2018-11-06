#ifndef LIMITERRIEMANNWENOS_H
#define LIMITERRIEMANNWENOS_H

#include "Limiter.h"

class LimiterRiemannWENOS : public Limiter
{

public:
    //- Construct by indicator and problem
    LimiterRiemannWENOS(const Indicator& ind, Problem& prb) : Limiter(ind, prb) {}

    //- Limit solution gradients
    virtual void limit(std::vector<numvector<double, 5 * nShapes> >& alpha) override;
};

#endif // LIMITERRIEMANNWENOS_H
