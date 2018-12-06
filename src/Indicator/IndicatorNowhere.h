#ifndef INDICATORNOWHERE_H
#define INDICATORNOWHERE_H

#include "Indicator.h"

class IndicatorNowhere : public Indicator
{
public:
    //- Constructor
    IndicatorNowhere (const Mesh& msh, const Physics& prb): Indicator(msh, prb) {}

    //- Check discontinuities
    virtual std::vector<int> checkDiscontinuities() const override;
};

#endif // INDICATORNOWHERE_H
