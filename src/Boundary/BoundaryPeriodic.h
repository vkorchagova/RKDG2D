#ifndef BOUNDARYPERIODIC_H
#define BOUNDARYPERIODIC_H

#include "Boundary.h"
#include <memory>

class Edge;

class BoundaryPeriodic : public Boundary
{
    std::shared_ptr<Edge> adjEdge;
    
public:
    BoundaryPeriodic(std::shared_ptr<Edge> e) : Boundary(), adjEdge(e) { type = "periodic";}
    ~BoundaryPeriodic() {}

    //- Apply boundary condition
    numvector<double, 5> applyBoundary(const numvector<double, 5>& solLeft = {0.0, 0.0, 0.0, 0.0, 0.0}, const Point& n = Point({0.0,0.0}), int numGP = 0) const override;
};

#endif // BOUNDARYPERIODIC_H
