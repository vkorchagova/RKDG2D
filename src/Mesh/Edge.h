/// ------------------------------
/// Edge class
/// ------------------------------
/// Consists of two nodes
/// Know its Gauss points
/// Can evaluate numerical fluxes through gauss points
/// ------------------------------

#ifndef EDGE_H
#define EDGE_H

#include "numvector.h"
#include "Cell.h"
#include "Point.h"

class Cell;

class Edge
{

public:


    /// geometric variables

    //- Two nodes define edge
    numvector<Point*, 2> nodes;

    //- Neighbour cells for edge
    numvector<Cell*, 2> neibCells;


    /// RKDG variables

    //- Number of gauss points for edge
    static const int nGP = 2;

    //- Local numerical fluxes for edge
    numvector<double,2> localFluxes;

public:

    //- Default constructor
    Edge();

    //- Copy constructor
    Edge(const Edge&);

    //- Overloaded "=" operator
    Edge& operator=(const Edge&);

    //- Destructor
    ~Edge();


    /// RKDG methods

    //- Calculate local fluxes for edge
    void getLocalFluxes();

}; // for Edge

#endif // EDGE_H
