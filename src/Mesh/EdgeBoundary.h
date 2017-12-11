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

class EdgeBoundary : public Edge
{
    /// geometric variables

    static const int nNeighbourCells = 1;

public:
    //- Default constructor
    EdgeBoundary() : Edge() {}

    //- Construct using two nodes
    EdgeBoundary(Point* p1, Point* p2) : Edge(p1, p2) { neibCells.reserve(nNeighbourCells); }

    //- Copy constructor
    EdgeBoundary(const EdgeBoundary& e) : Edge(e) {}

    //- Destructor
    ~EdgeBoundary() {}

    //- Apply boundary conditions
    virtual numvector<double, 5> applyBoundary(numvector<double, 5>& sol);
};

#endif // EDGEBOUNDARY_H
