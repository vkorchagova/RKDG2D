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
    /// Linear weight for WENO_S algorithm
    const double g = 0.001;

public:

    /// Constructor
    LimiterWENOS(
        const std::vector<std::shared_ptr<Cell>>& cells, 
        const Solution& sln,
        const Physics& phs);

    /// Destruct
    ~LimiterWENOS() {};

    /// Limit solution gradients
    virtual void limit(std::vector<numvector<double, dimS> >& SOL) override;
};

#endif // LIMITERWENOS_H
