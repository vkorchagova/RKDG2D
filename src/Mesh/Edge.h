/// ------------------------------
/// Edge
/// ------------------------------
/// Abstract class
///
/// Parameters:
/// -   nodes;
/// -   Gauss points
/// -   local fluxes in Gauss poins
///
/// Methods:
/// -   get local fluxes (public, pure virtual)
/// ------------------------------

#ifndef EDGE_H
#define EDGE_H

#include "numvector.h"
#include "Point.h"
#include "Cell.h"
#include "Flux.h"
#include <functional>
#include <iostream>

class Edge
{
private:

public:
    //- Number of gauss points for edge
    static const int nGP = 2;

    //- Gauss points
    numvector<Point, nGP> gPoints;

    //- Weights for integration
    numvector<double, nGP> gWeights;

    //- Jacobian
    double J;

    /// geometric variables

    //- Two nodes define edge
    numvector<std::shared_ptr<Point>, 2> nodes;

    //- Neighbour cells for edge: 2 for internal, 1 for boundary
    std::vector<std::shared_ptr<Cell>> neibCells;

    //- Normal to edge
    Point n;

    /// RKDG variables

    //- Local numerical fluxes for edge
    numvector<numvector<double, 5>, nGP> localFluxes;

public:

    //- Default constructor
    Edge() {}

    //- Construct using two nodes
    Edge(const Point& p1, const Point& p2);

    //- Copy constructor
    //Edge (const Edge&) = delete;

    //- Overloaded "=" operator
    //Edge& operator=(const Edge&) = delete;

    //- Destructor
    virtual ~Edge() = default;

    /// RKDG methods
    
    //- Calculate local fluxes for edge
    virtual void getLocalFluxesHor(const Flux& flux) = 0;
    virtual void getLocalFluxesVer(const Flux& flux) = 0;

    virtual void setBoundaryFunction(const numvector<double, 5>& bc) = 0;

    //- Calculate 1D integral through edge
    numvector<double, 5 * nShapes> boundaryIntegral(const std::shared_ptr<Cell>& cell) const;

}; // for Edge

#endif // EDGE_H
