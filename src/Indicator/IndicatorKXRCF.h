#ifndef INDICATORKXRCF_H
#define INDICATORKXRCF_H

#include "Indicator.h"
#include <vector>


class IndicatorKXRCF : public Indicator
{

public:

    //- Constructor
    IndicatorKXRCF (const Mesh2D& msh, const Problem& prb): Indicator (msh,prb) {}
    
    
    //- Check discontinuities
    virtual std::vector<int> checkDiscontinuities() const override;
};



#endif
