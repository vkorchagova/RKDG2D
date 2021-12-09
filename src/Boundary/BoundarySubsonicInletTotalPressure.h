#ifndef BOUNDARYSUBSONICINLETTOTALPRESSURE_H
#define BOUNDARYSUBSONICINLETTOTALPRESSURE_H

#include "Boundary.h"

/// 
/// Boundary condition: constant value for all conservative variables
/// 

class BoundarySubsonicInletTotalPressure : public Boundary
{
    double pTot;
    double TTot;

public:

    /// Constructor
    BoundarySubsonicInletTotalPressure(const Patch& p, const Physics& _phs, double _pTot, double _TTot) : Boundary(p, _phs), pTot(_pTot), TTot(_TTot) { type = "constant"; }

    // Get solution in outer side
    virtual numvector<double, dimPh> getSolOuter(
        const numvector<double, dimPh>& solInner, 
        const Point& n = Point({0.0,0.0})
    ) const override;

};

#endif // BoundarySubsonicInletTotalPressure_H
