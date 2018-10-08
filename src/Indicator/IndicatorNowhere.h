#ifndef INDICATORNOWHERE_H
#define INDICATORNOWHERE_H

#include "Indicator.h"

///
/// Optimistic indicator
///
/// Sets no one cell is troubled
///

class IndicatorNowhere : public Indicator
{
public:
    /// Constructor
    IndicatorNowhere (const Mesh2D& msh, const Problem& prb): Indicator(msh, prb) {}

    /// Check discontinuities
    virtual std::vector<int> checkDiscontinuities() const override;
};

#endif // INDICATORNOWHERE_H
