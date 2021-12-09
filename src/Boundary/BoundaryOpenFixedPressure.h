#ifndef BOUNDARYOPENFIXEDPRESSURE_H
#define BOUNDARYOPENFIXEDPRESSURE_H

#include "Boundary.h"

/// 
/// Open boundary condition
/// 

class BoundaryOpenFixedPressure : public Boundary
{
    double pressure;

public:

    /// Default constructor
    BoundaryOpenFixedPressure(const Patch& p, const Physics& _phs, double _pres = 101325) : Boundary(p, _phs), pressure(_pres) { type = "open_fixed_pressure";}

    // Get solution in outer side
    virtual numvector<double, dimPh> getSolOuter(
        const numvector<double, dimPh>& solInner, 
        const Point& n = Point({0.0,0.0})
    ) const override;

};

#endif // BoundaryOpenFixedPressure_H
