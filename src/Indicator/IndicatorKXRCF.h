#ifndef INDICATORKXRCF_H
#define INDICATORKXRCF_H

#include "Indicator.h"
#include <vector>


class IndicatorKXRCF : public Indicator
{

public:

    //- Constructor
    IndicatorKXRCF (const Mesh2D& msh): Indicator (msh) {};
    
    
    //- Check discontinuities
    virtual std::vector<double> checkDiscontinuities() const override;
};



#endif
