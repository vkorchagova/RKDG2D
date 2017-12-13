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
    EdgeInternal(const Flux& flux_) : Edge(flux_) {}

    //- Construct using two nodes
    EdgeInternal(const Point& p1, const Point& p2, const Flux& flux_) : Edge(p1, p2, flux_) { neibCells.reserve(nNeighbourCells); }

    //- Destructor
    ~EdgeInternal() {}


    virtual void getLocalFluxes() const override;
};

#endif // EDGEINTERNAL_H
