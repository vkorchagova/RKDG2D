#ifndef INDICATOREVERYWHERE_H
#define INDICATOREVERYWHERE_H

#include "Indicator.h"

class IndicatorEverywhere : public Indicator
{
public:
    //- Constructor
    IndicatorEverywhere (const Mesh& msh, const Physics& prb): Indicator (msh, prb) {}

    //- Check discontinuities
    virtual std::vector<int> checkDiscontinuities() const override;
};

#endif // INDICATOREVERYWHERE_H
