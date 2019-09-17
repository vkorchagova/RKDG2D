#ifndef EDGE_H
#define EDGE_H

#include "numvector.h"
#include "Cell.h"
#include <iostream>
#include <memory>

class Cell;

///
/// Edge 
///

class Edge
{
private:

    /// Lenght of edge ??TODO PERIMETER IN CELL
    double length;
  

public:

    /// Edge ID
    int number;

    /// Number of gauss points for edge
    int nGP;

    /// Gauss points
    std::vector<Point> gPoints;

    /// Weights for integration
    std::vector<double> gWeights;

    /// Jacobian
    double J;

    /// Two nodes define edge
    std::vector<std::shared_ptr<Point>> nodes;

    /// Normal to edge
    Point n;

    /// Neighbour cells for edge
    std::vector<std::shared_ptr<Cell>> neibCells;

    /// Copy constructor
    Edge(const Edge&) = default;

    /// Construct using two nodes
    Edge(const std::vector<std::shared_ptr<Point>> &p);

    /// Destructor
    virtual ~Edge() = default;

    /// Get length
    double getLength() const { return length; }

}; // for Edge

#endif // EDGE_H
