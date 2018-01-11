#ifndef EDGEBOUNDARYOPEN_H
#define EDGEBOUNDARYOPEN_H

#include "EdgeBoundary.h"


class EdgeBoundaryOpen : public EdgeBoundary
{
public:
    //- Default constructor
    EdgeBoundaryOpen() : EdgeBoundary()  {}

    //- Construct using two nodes
    EdgeBoundaryOpen(const Point& p1, const Point& p2) : EdgeBoundary(p1, p2) {}

    //- Destructor
    virtual ~EdgeBoundaryOpen() = default;

    //// RKDG methods

    //- Set BC function
    virtual void setBoundaryFunction(const numvector<double, 5>& bc) override { }

    //- Apply boundary condition
    numvector<double, 5> applyBoundary(const numvector<double, 5>& solLeft = {0.0, 0.0, 0.0, 0.0, 0.0}) const override;

};

#endif // EDGEBOUNDARYOPEN_H
