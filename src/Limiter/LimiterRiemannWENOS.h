#ifndef LIMITERRIEMANNWENOS_H
#define LIMITERRIEMANNWENOS_H

#include "LimiterWENOS.h"

///
/// WENO_S limiter applied to Riemann invariants component-by-component
///

class LimiterRiemannWENOS : public LimiterWENOS
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

protected:
    
    /// Limit solution gradients
    virtual numvector<double, dimS> limitation(const std::vector<std::shared_ptr<Cell>>& stencil) override;

    /// Choose stencil for defined cell
    virtual std::vector<std::shared_ptr<Cell>> getStencilFor(const std::shared_ptr<Cell>& cell) override;

public:

    /// Constructor
    LimiterRiemannWENOS(
        const Mesh& msh,
        Solution& sln,
        const Physics& phs,
        const Indicator& ind) : LimiterWENOS(msh, sln, phs, ind) {}

    /// Destructor
    ~LimiterRiemannWENOS() {};

};

#endif // LIMITERRIEMANNWENOS_H
