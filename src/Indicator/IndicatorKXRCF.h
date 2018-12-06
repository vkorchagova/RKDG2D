#ifndef INDICATORKXRCF_H
#define INDICATORKXRCF_H

#include "Indicator.h"

class IndicatorKXRCF : public Indicator
{

    //- Get mass fluxes
    numvector<double, 4> massFlux(const Edge& edge, const Cell& cell) const;

public:

    //- Constructor
    IndicatorKXRCF (const Mesh& msh, const Physics& prb): Indicator (msh, prb) {}
    
    
    //- Check discontinuities
    virtual std::vector<int> checkDiscontinuities() const override;


};



#endif
