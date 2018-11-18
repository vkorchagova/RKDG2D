#ifndef EDGE_H
#define EDGE_H

#include "numvector.h"
#include "defs.h"
//#include "Point.h"
#include "Cell.h"
//#include <functional>
#include <iostream>
#include <memory>

class Cell;

class Edge
{
public:

    /// geometric variables

    //- Number of gauss points for edge
    int nGP;

    //- Gauss points
    std::vector<Point> gPoints;

    //- Weights for integration
    std::vector<double> gWeights;

    //- Lenght of edge
    double length;

    //- Jacobian
    double J;


    //- Two nodes define edge
    //std::vector<std::unique_ptr<Point>> nodes;
    std::vector<std::reference_wrapper<Point>> nodes;

    //- Neighbour cells for edge
    //std::vector<std::unique_ptr<Cell>> neibCells;
    std::vector<std::reference_wrapper<Cell>> neibCells;

    //- Normal to edge
    Point n;

    //- Default constructor
    //Edge() {}
    
    //- Copy constructor
    //Edge(const Edge&) = delete;

    //- Construct using two nodes
    //Edge(const Point& p1, const Point& p2);
    Edge(const std::vector<std::reference_wrapper<Point>> &p);

    //- Destructor
    virtual ~Edge() = default;

    //- Get length
    double getLength() const { return length; }

}; // for Edge

#endif // EDGE_H
