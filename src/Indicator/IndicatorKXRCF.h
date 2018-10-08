#ifndef INDICATORKXRCF_H
#define INDICATORKXRCF_H

#include "Indicator.h"

///
/// KXRCF indicator
///
/// Check troubled cells by algorithm of Krivodonova et al.
///

class IndicatorKXRCF : public Indicator
{

    /// Get mass fluxes through each edge
    numvector<double, 4> massFlux(const Edge& edge, const Cell& cell) const;

public:

    /// Constructor
    IndicatorKXRCF (const Mesh2D& msh, const Problem& prb): Indicator (msh, prb) {}
    
    /// Check discontinuities
    virtual std::vector<int> checkDiscontinuities() const override;
};

#endif // INDICATORKXRCF_H
