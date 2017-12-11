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
    numvector<double, 5> infty;

public:
    //- Default constructor
    EdgeBoundaryInfty() : EdgeBoundary() {}

    //- Construct using two nodes
    EdgeBoundaryInfty(Point* p1, Point* p2) : EdgeBoundary(p1, p2) {}

    //- Copy constructor
    EdgeBoundaryInfty(const EdgeBoundaryInfty& e) : EdgeBoundary(e) {}

    //- Destructor
    ~EdgeBoundaryInfty() {}

    //// RKDG methods

    numvector<double, 5> applyBoundary();

    void setCondition(numvector<double, 5>& infValues) { infty = infValues; }
};

#endif // EDGEBOUNDARYINFTY_H
