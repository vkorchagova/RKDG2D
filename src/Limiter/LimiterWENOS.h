#ifndef LIMITERWENOS_H
#define LIMITERWENOS_H

#include "Limiter.h"

///
/// WENO_S limiter applied to conservative variables component-by-component
///

/// Log file to save data
extern std::ofstream logger;

class LimiterWENOS : public Limiter
{

protected:

    /// Linear weight for WENO_S algorithm
    const double g = 0.001;
    
    /// Limit solution gradients
    virtual numvector<double, dimS> limitation(const std::vector<std::shared_ptr<Cell>>& stencil) override;

    /// Choose stencil for defined cell
    virtual std::vector<std::shared_ptr<Cell>> getStencilFor(const std::shared_ptr<Cell>& cell) override;

    /// Limit polynoms
    numvector<double, dimS> limitP(
        const std::vector<std::shared_ptr<Cell>>& stencil, 
        const std::vector<numvector<double, dimS>>& p);

public:

    /// Constructor
    LimiterWENOS(
        const Mesh& msh,
        Solution& sln,
        const Physics& phs,
        const Indicator& ind,
        Buffers& buf);

    /// Destructor
    ~LimiterWENOS() {};
};

#endif // LIMITERWENOS_H
