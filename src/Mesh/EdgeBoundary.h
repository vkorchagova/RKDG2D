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
    EdgeBoundary(const Flux& flux_) : Edge(flux_) {};

    //- Construct using two nodesD
    EdgeBoundary(const Point& p1, const Point& p2, const Flux& flux_) : Edge(p1, p2, flux_) { neibCells.reserve(nNeighbourCells); }

    //- Destructor
    virtual ~EdgeBoundary() {}

    //- Apply boundary conditions
    virtual numvector<double, 5> applyBoundary(const numvector<double, 5>& sol = {0.0, 0.0, 0.0, 0.0, 0.0}) const = 0;

    virtual void getLocalFluxes() const override;
};

#endif // EDGEBOUNDARY_H
