/// ------------------------------
/// EdgeInternal
/// ------------------------------
/// Base class: Edge
///
/// Parameters:
/// -   two neighbour cells
///
/// Methods:
/// -   get local fluxes (public)
/// ------------------------------

#ifndef EDGEINTERNAL_H
#define EDGEINTERNAL_H

#include "Edge.h"

class EdgeInternal : public Edge
{

    /// geometric variables

    static const int nNeighbourCells = 2;


public:

    //- Default constructor
    EdgeInternal() : Edge() {}

    //- Construct using two nodes
    EdgeInternal(const Point& p1, const Point& p2) : Edge(p1, p2) { neibCells.reserve(nNeighbourCells); }

    //- Destructor
    virtual ~EdgeInternal() = default;

    virtual void setBoundaryFunction(const numvector<double, 5>& bc) override {}

    //- Calculate local fluxes in gauss points
    virtual void getLocalFluxes(const Flux& flux) override;
};

#endif // EDGEINTERNAL_H
