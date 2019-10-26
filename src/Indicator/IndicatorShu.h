#ifndef INDICATORSHU_H
#define INDICATORSHU_H

#include "Indicator.h"

///
/// Indicator SHU
///

class IndicatorShu : public Indicator
{
    // Indication constant
    const double Ck = 0.03;

    // Unzero-denom
    double eps = 1e-6;

public:

    /// Constructor
    IndicatorShu(const Mesh& msh, const Solution& sln) : Indicator(msh, sln) {}

    /// Destructor
    ~IndicatorShu() {}

    /// Check discontinuities
    virtual std::vector<int> checkDiscontinuities() const override;
};

#endif // INDICATORSHU_H
