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
    EdgeInternal(Point* p1, Point* p2) : Edge(p1, p2) { neibCells.reserve(nNeighbourCells); }

    //- Copy constructor
    EdgeInternal(const EdgeInternal& e) : Edge(e) {}

    //- Destructor
    ~EdgeInternal() {}


    void getLocalFluxes() {}
};

#endif // EDGEINTERNAL_H
