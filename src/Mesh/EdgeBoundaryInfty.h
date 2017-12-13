/// ------------------------------
/// EdgeBoundaryInfty
/// ------------------------------
/// Base class: EdgeBoundary
///
/// ABoundary condition on infinity
///
/// Parameters:
/// -   function for infinity state
///
/// Methods:
/// -   get local fluxes (public)
/// -   apply boundary condition (private)
/// ------------------------------



#ifndef EDGEBOUNDARYINFTY_H
#define EDGEBOUNDARYINFTY_H

#include "EdgeBoundary.h"


class EdgeBoundaryInfty : public  EdgeBoundary
{

public:

    numvector<double, 5> infty;

public:
    //- Default constructor
    EdgeBoundaryInfty(const Flux& flux_) : EdgeBoundary(flux_) {};

    //- Construct using two nodes
    EdgeBoundaryInfty(const Point& p1, const Point& p2, const Flux& flux_) : EdgeBoundary(p1, p2, flux_) {};

    //- Destructor
    ~EdgeBoundaryInfty() {};

    //// RKDG methods

    numvector<double, 5> applyBoundary(const numvector<double, 5>& solLeft = {0.0, 0.0, 0.0, 0.0, 0.0}) const override;

    //void setCondition(numvector<double, 5>& infValues) { infty = infValues; }
};

#endif // EDGEBOUNDARYINFTY_H
