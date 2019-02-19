#ifndef EDGE_H
#define EDGE_H

#include "numvector.h"
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


    int number;

    /// Number of gauss points for edge
    int nGP;

    /// Gauss points
    std::vector<Point> gPoints;

    /// Weights for integration
    std::vector<double> gWeights;

    /// Lenght of edge
    double length;

    /// Jacobian
    double J;


    /// Two nodes define edge
    //std::vector<std::shared_ptr<Point>> nodes;
    std::vector<std::shared_ptr<Point>> nodes;

    /// Neighbour cells for edge
    //std::vector<std::shared_ptr<Cell>> neibCells;
    std::vector<std::shared_ptr<Cell>> neibCells;

    /// Normal to edge
    Point n;

    /// Default constructor
    //Edge() {}
    
    /// Copy constructor
    Edge(const Edge&) = default;

    /// Construct using two nodes
    //Edge(const Point& p1, const Point& p2);
    Edge(const std::vector<std::shared_ptr<Point>> &p);

    /// Destructor
    virtual ~Edge() = default;

    /// Get length
    double getLength() const { return length; }

    /// Check edge is equal to given
    bool isEqual(const Edge& e) const;
    bool isEqual(const std::shared_ptr<Edge>& e) const { return isEqual(*e); }

}; // for Edge

#endif // EDGE_H
