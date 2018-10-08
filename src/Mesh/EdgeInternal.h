#ifndef EDGEINTERNAL_H
#define EDGEINTERNAL_H

#include "Edge.h"

///
/// EdgeInternal class
///
/// Parameters:
/// -   two neighbour cells
///
/// Methods:
/// -   get local fluxes (public)

class EdgeInternal : public Edge
{
    /// number of neighbour cells
    static const int nNeighbourCells = 2;


public:

    /// Default constructor
    EdgeInternal() : Edge() {}

    /// Construct using two nodes
    EdgeInternal(const Node& p1, const Node& p2) : Edge(p1, p2) { neibCells.reserve(nNeighbourCells); }

    /// Destructor
    virtual ~EdgeInternal() = default;

    /// Set type of boundary condition - just empty stopper - trace of EdgeBoundaryDiagProjection...
    virtual void setBoundary(const std::shared_ptr<Boundary>& bound) override {}

    /// Calculate local fluxes in gauss points
    virtual void getLocalFluxes(const Flux& flux) override;

    /// Compute max speed on edge
    virtual void getMaxUL() override;
};

#endif // EDGEINTERNAL_H
