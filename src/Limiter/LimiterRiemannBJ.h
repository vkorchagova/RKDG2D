#ifndef LIMITERRIEMANNBJ_H
#define LIMITERRIEMANNBJ_H

#include "LimiterBJ.h"

/// Log file to save data
extern std::ofstream logger;

///
/// BJ limiter applied to conservative variables component-by-component in gauss points of edges
///

class LimiterRiemannBJ : public LimiterBJ
{
                       
    /// Turn to conservative variables from Riemann
    numvector<double, dimS> riemannToConservative
    (
        const numvector<double, dimS>& alpha
    ) const;

    /// Turn to Riemann variables from conservative
    numvector<double, dimS> conservativeToRiemann
    (
        const numvector<double, dimS>& alpha
    ) const;

    numvector<numvector<double, dimPh>, dimPh>  L;
    numvector<numvector<double, dimPh>, dimPh>  R;

    numvector<double, dimS> limitChar(const std::vector<std::shared_ptr<Cell>>& stencil, const std::vector<numvector<double, dimS>>& pInv);

protected:

    /// Limit solution gradients
    ///
    /// @param alpha vector of solution coeffitients in all cells which should be limited
    virtual numvector<double, dimS> limitation(const std::vector<std::shared_ptr<Cell>>& stencil) override;

public:

    /// Constructor 
    LimiterRiemannBJ(
        const Mesh& mesh, 
        Solution& sln,
        const Physics& phs,
        const Indicator& ind,
        Buffers& _buf) : LimiterBJ(mesh, sln, phs, ind, _buf) {}; 

    /// Destructor
    ~LimiterRiemannBJ() {};   
};

#endif // LIMITERRIEMANNBJ_H
