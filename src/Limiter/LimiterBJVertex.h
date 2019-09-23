#ifndef LIMITERBJVERTEX_H
#define LIMITERBJVERTEX_H

#include "Limiter.h"

///
/// BJ limiter applied to conservative variables component-by-component by gauss points of edges and cell tertices additionally
///

class LimiterBJVertex : public Limiter
{
    /// Get AlphaL
    numvector<double, dimPh> getYMin(
                                const std::shared_ptr<Cell>& cell,
                                const numvector<double, dimPh>& mI,
                                const numvector<double, dimPh>& MI,
                                const numvector<double, dimPh>& uMean
                            );

protected:

    /// Limit solution gradients
    ///
    /// @param alpha    vector of solution coeffitients in all cells which should be limited
    virtual numvector<double, dimS> limitation(const std::vector<std::shared_ptr<Cell>>& stencil) override;

    /// Choose stencil for defined cell
    virtual std::vector<std::shared_ptr<Cell>> getStencilFor(const std::shared_ptr<Cell>& cell) override;

public:
    /// Construct
    LimiterBJVertex(
        const Mesh& msh,
        Solution& sln,
        const Physics& phs,
        const Indicator& ind);

    /// Destructor
    ~LimiterBJVertex() {};
};

#endif // LIMITERBJVERTEX_H
