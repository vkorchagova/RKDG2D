/// ------------------------------
/// EdgeBoundary
/// ------------------------------
/// Base class: Edge
///
/// Abstract class
///
/// Parameters:
/// -   one neighbour cell
///
/// Methods:
/// -   get local fluxes (public)
/// -   apply boundary condition (private, pure virtual)
/// ------------------------------


#ifndef EDGEBOUNDARY_H
#define EDGEBOUNDARY_H

#include "Edge.h"
#include "Boundary.h"

class EdgeBoundary : public Edge
{
    /// geometric variables

    static const int nNeighbourCells = 1;

public:

    //- Boundary condition
    std::shared_ptr<Boundary> bc;

    //- Default constructor
    EdgeBoundary() : Edge() {}

    //- Construct using two nodes
    EdgeBoundary(const Point& p1, const Point& p2) : Edge(p1, p2) { neibCells.reserve(nNeighbourCells); }

    //- Destructor
    virtual ~EdgeBoundary() = default;

    // virtual void setMeshPointer(const Mesh2D& msh) {};

    //- Set boundary condition
    void setBoundary(const std::shared_ptr<Boundary>& bound) {bc = bound;}

    //- Calculate local fluxes in gauss points
    virtual void getLocalFluxes(const Flux& flux) override;

    //- Compute max speed on edge
    virtual void getMaxUL() override;
};

#endif // EDGEBOUNDARY_H
