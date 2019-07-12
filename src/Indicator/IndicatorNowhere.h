#ifndef INDICATORNOWHERE_H
#define INDICATORNOWHERE_H

#include "Indicator.h"

///
/// Pure optimistic indicator of problem cells: all cells considered as good
///

class IndicatorNowhere : public Indicator
{
public:

    /// Constructor
    IndicatorNowhere(const Mesh& msh) : Indicator(msh) {}

    /// Check discontinuities
    virtual std::vector<int> checkDiscontinuities() const override;
};

#endif // INDICATORNOWHERE_H
