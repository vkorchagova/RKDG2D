#ifndef BOUNDARYOPENTOTALPRESSURE_H
#define BOUNDARYOPENTOTALPRESSURE_H

#include "Boundary.h"

/// 
/// Open boundary condition
/// 

class BoundaryOpenTotalPressure : public Boundary
{
    double pTotal;
    double TFix;

public:

    /// Default constructor
    BoundaryOpenTotalPressure(const Patch& p, const Physics& _phs, double _pres = 101325, double _TFix = 293.15) : Boundary(p, _phs), pTotal(_pres), TFix(_TFix) { type = "open_total_pressure";}

    // Get solution in outer side
    virtual numvector<double, dimPh> getSolOuter(
        const numvector<double, dimPh>& solInner, 
        const Point& n = Point({0.0,0.0})
    ) const override;

};

#endif // BoundaryOpenTotalPressure_H
