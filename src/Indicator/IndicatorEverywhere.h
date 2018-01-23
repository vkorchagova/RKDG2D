#ifndef INDICATOREVERYWHERE_H
#define INDICATOREVERYWHERE_H

#include "Indicator.h"

class IndicatorEverywhere : public Indicator
{
public:
    //- Constructor
    IndicatorEverywhere (const Mesh2D& msh): Indicator (msh) {}

    //- Check discontinuities
    virtual std::vector<int> checkDiscontinuities() const override;
};

#endif // INDICATOREVERYWHERE_H
