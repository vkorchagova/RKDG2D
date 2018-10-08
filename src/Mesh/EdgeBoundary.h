#ifndef EDGEBOUNDARY_H
#define EDGEBOUNDARY_H

#include "Edge.h"
#include "Boundary.h"

/// Abstract EdgeBoundary class
///
/// Parameters:
/// -   one neighbour cell
///
/// Methods:
/// -   get local fluxes (public)

class EdgeBoundary : public Edge
{
    //- geometric variables

    /// number of neighbour cells
    static const int nNeighbourCells = 1;

public:

    /// Boundary condition
    std::shared_ptr<Boundary> bc;

    /// Default constructor
    EdgeBoundary() : Edge() {}

    /// Construct using two nodes
    EdgeBoundary(const Node& p1, const Node& p2) : Edge(p1, p2) { neibCells.reserve(nNeighbourCells); }

    /// Destructor
    virtual ~EdgeBoundary() = default;

    /// Set boundary condition
    void setBoundary(const std::shared_ptr<Boundary>& bound) {bc = bound;}

    /// Calculate local fluxes in gauss points
    virtual void getLocalFluxes(const Flux& flux) override;

    /// Compute max speed on edge
    virtual void getMaxUL() override;
};

#endif // EDGEBOUNDARY_H
