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

    /// geometric variables

    //- Two nodes define edge
    numvector<const Point*, 2> nodes;

    //- Neighbour cells for edge: 2 for internal, 1 for boundary
    std::vector<Cell*> neibCells;

    const Flux& flux;

    /// RKDG variables

    //- Local numerical fluxes for edge
    numvector<numvector<double, 5>, nGP> localFluxes;


public:

    //- Default constructor
    Edge(const Flux& flux_) : flux(flux_) {}

    //- Construct using two nodes
    Edge(const Point& p1, const Point& p2, const Flux& flux_);

    //- Destructor
    virtual ~Edge() {}

    /// RKDG methods
    
    //- Calculate local fluxes for edge
    virtual void getLocalFluxes() const = 0;
    
    //- Set flux
    void setFlux (const Flux& flux_) const { flux = flux_; };

    //- Calculate 1D integral through edge
    numvector<double, 5> boundaryIntegral(const std::function<double(const Point&)>& phi) const;

}; // for Edge

#endif // EDGE_H
