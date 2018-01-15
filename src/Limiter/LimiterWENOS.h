#ifndef LIMITERWENOS_H
#define LIMITERWENOS_H

#include "Limiter.h"


class LimiterWENOS : public Limiter
{
    void limitSep(std::vector<numvector<double, 5 * nShapes> >& alpha);
    void limitTog(std::vector<numvector<double, 5 * nShapes> >& alpha);
    void limitX(std::vector<numvector<double, 5 * nShapes> >& alpha);

public:
    //- Construct by indicator and problem
    LimiterWENOS(const Indicator& ind, Problem& prb) : Limiter(ind, prb) {}

    //- Limit solution gradients
    virtual void limit(std::vector<numvector<double, 5 * nShapes> >& alpha) override;
};

#endif // LIMITERWENOS_H
