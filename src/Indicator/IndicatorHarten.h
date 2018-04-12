#ifndef INDICATORHARTEN_H
#define INDICATORHARTEN_H

#include "Indicator.h"

class IndicatorHarten : public Indicator
{


public:

    //- Constructor
    IndicatorHarten (const Mesh2D& msh, const Problem& prb): Indicator (msh, prb) {}

    //- Check discontinuities
    virtual std::vector<int> checkDiscontinuities() const override;
};

#endif // INDICATORHARTEN_H
