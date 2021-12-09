#ifndef BOUNDARYSUBSONICINLETFIXEDPRESSURE_H
#define BOUNDARYSUBSONICINLETFIXEDPRESSURE_H

#include "Boundary.h"

/// 
/// Boundary condition: constant value for all conservative variables
/// 

class BoundarySubsonicInletFixedPressure : public Boundary
{
    double pFix;
    double TTot;

public:

    /// Constructor
    BoundarySubsonicInletFixedPressure(const Patch& p, const Physics& _phs, double _pFix, double _TTot) : Boundary(p, _phs), pFix(_pFix), TTot(_TTot) { type = "constant"; }

    // Get solution in outer side
    virtual numvector<double, dimPh> getSolOuter(
        const numvector<double, dimPh>& solInner, 
        const Point& n = Point({0.0,0.0})
    ) const override;

};

#endif // BoundarySubsonicInletFixedPressure_H
