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
    EdgeBoundaryInfty() : EdgeBoundary()  {};

    //- Construct using two nodes
    EdgeBoundaryInfty(const Point& p1, const Point& p2) : EdgeBoundary(p1, p2) {};

    //- Destructor
    virtual ~EdgeBoundaryInfty() = default;

    //// RKDG methods

    //- Set BC function
    virtual void setBoundaryFunction(const numvector<double, 5>& bc) override { infty = bc; }

    numvector<double, 5> applyBoundary(const numvector<double, 5>& solLeft = {0.0, 0.0, 0.0, 0.0, 0.0}) const override;
};

#endif // EDGEBOUNDARYINFTY_H
