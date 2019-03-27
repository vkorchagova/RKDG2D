#ifndef LIMITERWENOS_H
#define LIMITERWENOS_H

#include "Limiter.h"

/// Log file to save data
extern std::ofstream logger;

class LimiterWENOS : public Limiter
{
    const double g = 0.001;

public:
    //- Construct by indicator and problem
    //LimiterWENOS(const Indicator& ind, Physics& phys, Solution& sol, Mesh& mesh) : Limiter(ind, phys, sol, mesh) {}
    LimiterWENOS(
        const std::vector<std::shared_ptr<Cell>>& cells, 
        const Solution& sln,
        const Physics& phs);

    ~LimiterWENOS() {};

    //- Limit solution gradients
    virtual void limit(std::vector<numvector<double, dimS> >& SOL) override;
};

#endif // LIMITERWENOS_H
