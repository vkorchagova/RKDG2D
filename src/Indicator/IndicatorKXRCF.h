#ifndef INDICATORKXRCF_H
#define INDICATORKXRCF_H

#include "Indicator.h"
#include <vector>


class IndicatorKXRCF : public Indicator
{

public:

    //- Constructor
    IndicatorKXRCF (const Mesh2D& msh, const std::vector<numvector<double, 5 * nShapes>>& coeffs): Indicator (msh, coeffs) {};
    
    
    //- Check discontinuities
    virtual std::vector<double> checkDiscontinuities() const override;
};



#endif
