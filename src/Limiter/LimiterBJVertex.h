#ifndef LIMITERBJVERTEX_H
#define LIMITERBJVERTEX_H

#include "Limiter.h"

///
/// BJ-VERTEX limiter applied to conservative variables component-by-component
///

class LimiterBJVertex : public Limiter
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
    LimiterBJVertex(
        const std::vector<std::shared_ptr<Cell>>& cells,
        const Solution& sln,
        const Physics& phs) : Limiter(cells, sln, phs) {}

    /// Limit solution gradients
    ///
    /// @param alpha    vector of solution coeffitients in all cells which should be limited
    virtual void limit(std::vector<numvector<double, dimS> >& alpha) override;
};

#endif // LIMITERBJVERTEX_H
