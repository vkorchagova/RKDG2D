#ifndef EDGEBOUNDARYSLIP_H
#define EDGEBOUNDARYSLIP_H

#include "EdgeBoundary.h"

class EdgeBoundarySlip : public EdgeBoundary
{
public:
    //- Default constructor
    EdgeBoundarySlip() : EdgeBoundary()  {}

    //- Construct using two nodes
    EdgeBoundarySlip(const Point& p1, const Point& p2) : EdgeBoundary(p1, p2) {}

    //- Destructor
    virtual ~EdgeBoundarySlip() = default;

    //// RKDG methods

    //- Set BC function
    virtual void setBoundaryFunction(const numvector<double, 5>& bc) override { }

    //- Apply boundary condition
    numvector<double, 5> applyBoundary(const numvector<double, 5>& solLeft = {0.0, 0.0, 0.0, 0.0, 0.0}) const override;
};

#endif // EDGEBOUNDARYSLIP_H
