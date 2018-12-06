#ifndef LIMITERBJ_H
#define LIMITERBJ_H

#include "Limiter.h"

///
/// WENO\_S limiter applied to conservative variables component-by-component
///

class LimiterBJ : public Limiter
{

public:
    ///- Construct by indicator and problem
    LimiterBJ(const Indicator& ind, Physics& phys, Solution& sol, Mesh& mesh) : Limiter(ind, phys, sol, mesh) {}
	
    /// Limit solution gradients
    ///
    /// @param SOL    vector of solution coeffitients in all cells which should be limited
    virtual void limit(std::vector<numvector<double, dimS> >& SOL) override;
};

#endif // LIMITERBJ_H
