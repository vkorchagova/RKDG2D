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
//#include "Mesh2D.h"

//class Mesh2D;

class EdgeBoundary : public Edge
{
    /// geometric variables

    static const int nNeighbourCells = 1;

public:
    //- Default constructor
    EdgeBoundary() : Edge() {}

    //- Construct using two nodes
    EdgeBoundary(const Point& p1, const Point& p2) : Edge(p1, p2) { neibCells.reserve(nNeighbourCells); }

    //- Destructor
    virtual ~EdgeBoundary() = default;

    virtual void setBoundaryFunction(const numvector<double, 5>& bc) override = 0;

   // virtual void setMeshPointer(const Mesh2D& msh) {};

    //- Apply boundary conditions
    virtual numvector<double, 5> applyBoundary(const numvector<double, 5>& sol = {0.0, 0.0, 0.0, 0.0, 0.0}) const = 0;

    //- Calculate local fluxes in gauss points
    virtual void getLocalFluxes(const Flux& flux) override;
};

#endif // EDGEBOUNDARY_H
