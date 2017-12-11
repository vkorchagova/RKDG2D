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
    numvector<Point*, nGP> nodes;

    //- Neighbour cells for edge: 2 for internal, 1 for boundary
    std::vector<Cell*> neibCells;

    Flux* flux;


    /// RKDG variables



    //- Local numerical fluxes for edge
    numvector<numvector<double, 5>, nGP> localFluxes;


public:

    //- Default constructor
    Edge() {}

    //- Construct using two nodes
    Edge(Point*, Point*);

    //- Copy constructor
    Edge(const Edge&);

    //- Overloaded "=" operator
    Edge& operator=(const Edge&);

    //- Destructor
    ~Edge() {}


    /// RKDG methods

    //- Set flux
    void setFlux(Flux& flx);

    //- Calculate local fluxes for edge
    virtual void getLocalFluxes() { std::cout << "i'm in base Edge \n";}

    //- Calculate 1D integral through edge
    numvector<double, 5> boundaryIntegral(std::function<double(const Point&)>& phi);

}; // for Edge

#endif // EDGE_H
