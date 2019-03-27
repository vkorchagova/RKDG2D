#ifndef LIMITERRIEMANNWENOS_H
#define LIMITERRIEMANNWENOS_H

#include "Limiter.h"

class LimiterRiemannWENOS : public Limiter
{

    /// Turn to conservative variables from Riemann
    numvector<double, dimS> riemannToConservative
    (
        const numvector<double, dimS>& alpha, 
        const numvector<numvector<double, dimPh>, dimPh>& R
    ) const;

    /// Turn to Riemann variables from conservative
    numvector<double, dimS> conservativeToRiemann
    (
        const numvector<double, dimS>& alpha, 
        const numvector<numvector<double, dimPh>, dimPh>& L
    ) const;

    // Limit Riemann variables in WENO_S manner
    //numvector<double, dimS> limitP
    //(
    //    const std::vector<std::shared_ptr<Cell>>& stenc, 
    //    const std::vector<numvector<double, dimS>>& alpha
    //);

public:
    //- Construct by indicator and problem
    LimiterRiemannWENOS(
        const std::vector<std::shared_ptr<Cell>>& cells, 
        const Solution& sln,
        const Physics& phs) : Limiter(cells, sln, phs) {}

    //- Limit solution gradients
    virtual void limit(std::vector<numvector<double, dimS> >& SOL) override;
};

#endif // LIMITERRIEMANNWENOS_H
