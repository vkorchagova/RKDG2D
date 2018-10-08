#ifndef INDICATOREVERYWHERE_H
#define INDICATOREVERYWHERE_H

#include "Indicator.h"

///
/// Pessimistic indicator
///
/// Sets all cells to be troubled
///

class IndicatorEverywhere : public Indicator
{
public:
    
    /// Constructor
    IndicatorEverywhere (const Mesh2D& msh, const Problem& prb): Indicator (msh, prb) {}

    /// Check discontinuities
    virtual std::vector<int> checkDiscontinuities() const override;
};

#endif // INDICATOREVERYWHERE_H
