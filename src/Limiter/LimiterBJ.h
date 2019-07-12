#ifndef LIMITERBJ_H
#define LIMITERBJ_H

#include "Limiter.h"

/// Log file to save data
extern std::ofstream logger;

///
/// BJ limiter applied to conservative variables component-by-component in gauss points of edges
///

class LimiterBJ : public Limiter
{
    /// Get AlphaL
    numvector<double, dimPh> getAlphaL(
                                const std::shared_ptr<Cell>& cell, 
								const numvector<double, dimPh>& mI,
								const numvector<double, dimPh>& MI,
								const numvector<double, dimPh>& uMean
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
