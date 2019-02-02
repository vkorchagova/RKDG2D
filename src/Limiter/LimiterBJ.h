#ifndef LIMITERBJ_H
#define LIMITERBJ_H

#include "Limiter.h"

///
/// WENO\_S limiter applied to conservative variables component-by-component
///

class LimiterBJ : public Limiter
{
    //- Get AlphaL
    numvector<double, dimPh> getAlphaL(
                                std::shared_ptr<Cell>& cell, 
                                numvector<double, dimPh>& mI, 
                                numvector<double, dimPh>& MI, 
                                numvector<double, dimPh>& uMean
                            );

public:
    /// Construct 
    LimiterBJ(
        const std::vector<std::shared_ptr<Cell>>& cells, 
        const Solution& sln,
        const Physics& phs) : Limiter(cells, sln, phs) {}

    /// Limit solution gradients
    ///
    /// @param alpha    vector of solution coeffitients in all cells which should be limited
    virtual void limit(std::vector<numvector<double, dimS> >& alpha) override;
};

#endif // LIMITERBJ_H
