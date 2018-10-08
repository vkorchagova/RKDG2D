#ifndef EDGE_H
#define EDGE_H

#include "numvector.h"
#include "defs.h"
#include "Point.h"
#include "Cell.h"
#include "Flux.h"
#include "Boundary.h"
#include <functional>
#include <iostream>

///
/// Abstract Edge class
///
/// Parameters:
/// -   nodes;
/// -   Gauss points
/// -   local fluxes in Gauss poins
///
/// Methods:
/// -   get local fluxes (public, pure virtual)
///

class Edge
{
private:

    /// Lenght of edge
    double length;

public:

    /// Number of gauss points for edge
    int nGP;

    /// Gauss points
    std::vector<Point> gPoints;

    /// Weights for integration
    std::vector<double> gWeights;

    /// Jacobian
    double J;

    //- geometric variables

    /// Edge number
    int number;

    /// Two nodes define edge
    numvector<std::shared_ptr<Node>, 2> nodes;

    /// Neighbour cells for edge: 2 for internal, 1 for boundary
    std::vector<std::shared_ptr<Cell>> neibCells;

    /// Normal to edge
    Point n;

    //- RKDG variables

    /// Local numerical fluxes for edge
    std::vector<numvector<double, 5>> localFluxes;

    /// Max velocity multiplied with edge length
    double uMaxL;

public:

    /// Default constructor
    Edge() {}

    /// Construct using two nodes
    Edge(const Node& p1, const Node& p2);

    /// Destructor
    virtual ~Edge() = default;

    /// Get length
    double getLength() const { return length; }
    
    //- RKDG methods
    
    /// Calculate local fluxes for edge
    virtual void getLocalFluxes(const Flux& flux) = 0;

    /// Set boundary condition --- implemented only for boundary edges
    virtual void setBoundary(const std::shared_ptr<Boundary>& bound) = 0;

    /// Calculate 1D integral through edge
    numvector<double, 5 * nShapes> boundaryIntegral(const std::shared_ptr<Cell>& cell) const;

    /// Compute mass flux through edge
    double getMassFlux(const std::shared_ptr<Cell> &cell) const;

    /// Compute max speed on edge
    virtual void getMaxUL() = 0;

    /// Return maxUL
    double maxUL() const {return uMaxL;}

}; // for Edge

#endif // EDGE_H
