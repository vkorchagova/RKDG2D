#ifndef LIMITERWENOS_H
#define LIMITERWENOS_H

#include "Limiter.h"


class LimiterWENOS : public Limiter
{

public:
    //- Construct by indicator and problem
    LimiterWENOS(const Indicator& ind, Physics& phys, Solution& sol, Mesh& mesh) : Limiter(ind, phys, sol, mesh) {}

    //- Limit solution gradients
    virtual void limit(std::vector<numvector<double, dimS> >& SOL) override;
};

#endif // LIMITERWENOS_H
