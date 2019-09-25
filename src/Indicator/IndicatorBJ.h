#ifndef INDICATORBJ_H
#define INDICATORBJ_H

#include "Indicator.h"

///
/// Indicator based on BJ limiter
///

class IndicatorBJ : public Indicator
{

    numvector<double, dimPh> getYMin(const std::shared_ptr<Cell>& cell, const numvector<double, dimPh>& mI, const numvector<double, dimPh>& MI, const numvector<double, dimPh>& uMean) const;

public:

    /// Constructor
    IndicatorBJ(const Mesh& msh, const Solution& sln) : Indicator(msh, sln) {}

    /// Destructor
    ~IndicatorBJ() {};

    /// Check discontinuities
    virtual std::vector<int> checkDiscontinuities() const override;
};

#endif // INDICATORBJ_H
