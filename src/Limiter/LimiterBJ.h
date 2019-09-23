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

private:

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
    /// @param alpha vector of solution coeffitients in all cells which should be limited
    virtual numvector<double, dimS> limitation(const std::vector<std::shared_ptr<Cell>>& stencil) override;

    /// Choose stencil for defined cell
    virtual std::vector<std::shared_ptr<Cell>> getStencilFor(const std::shared_ptr<Cell>& cell) override;

public:

    /// Constructor 
    LimiterBJ(
        const Mesh& mesh, 
        Solution& sln,
        const Physics& phs,
        const Indicator& ind); 

    /// Destructor
    ~LimiterBJ() {};   
};

#endif // LIMITERBJ_H
