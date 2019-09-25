#ifndef INDICATOREVERYWHERE_H
#define INDICATOREVERYWHERE_H

#include "Indicator.h"

///
/// Pure pessimistic indicator of problem cells: all cells considered as troubled
///

class IndicatorEverywhere : public Indicator
{
public:

    /// Constructor
    IndicatorEverywhere(const Mesh& msh, const Solution& sln) : Indicator(msh, sln) {}

    /// Destructor
    ~IndicatorEverywhere() {};

    /// Check discontinuities
    virtual std::vector<int> checkDiscontinuities() const override;
};

#endif // INDICATOREVERYWHERE_H
