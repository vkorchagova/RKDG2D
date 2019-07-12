#ifndef CELL_H
#define CELL_H

#include "numvector.h"

#include <functional>
#include <memory>
#include <math.h>
#include <algorithm>

#include "Point.h"
#include "Edge.h"

class Edge;

///
/// Сell 
///

class Cell
{

public:

    /// Cell ID
    int number;
    
    /// Number of Gauss points
    int nGP;

    /// Gauss points
    std::vector<Point> gPoints2D;
    
    /// Gauss weights
    std::vector<double> gWeights2D;
    
    /// Jacobian
    std::vector<double> J;
    
    /// Area of cell
    double area;
    
    /// Mass center of cell
    Point center;

    /// Gramian matrix
    std::vector<std::vector<double>> gramian;

    /// local [-1,1]x[-1,1] to global rectangular cell
    Point localToGlobal(const Point& localPoint) const;

    /// List of nodes in cell
    std::vector<std::shared_ptr<Point>> nodes;

    /// List of edges in cell
    std::vector<std::shared_ptr<Edge>> edges;

    /// Number of entities (nodes or edges)
    int nEntities;

    /// Neighbour cells
    std::vector<std::shared_ptr<Cell>> neibCells;

    /// Neighbour cells vertex
    std::vector<std::shared_ptr<Cell>> neibCellsVertex;

    /// Return area of cell
    double getArea() const { return area; }

    /// Return center of cell
    const Point& getCellCenter() const { return center; }

    /// Calculate element area
    void setArea();

    /// Set cell center
    void setCellCenter();

    /// Define Jacobian function
    void setJacobian();

    /// Set Gauss points
    void setGaussPoints();

    /// Construct cell using vectors of nodes and edges
    Cell(const std::vector<std::shared_ptr<Point>> &nodes, const std::vector<std::shared_ptr<Edge>> &edges);

    /// Destructor
    ~Cell() {}

    /// Calculate coordinates of cell nodes
    std::vector<Point> getCellCoordinates() const;

    /// Find neighbour cells
    void findNeighbourCells();

    /// check for a cell among vertex neighbor
    bool inNeibVertList(const std::shared_ptr<Cell>& cell, std::vector<std::shared_ptr<Cell>>& neibVC);

    /// check for common node
    bool hasCommonNode(const std::shared_ptr<Cell>& c1);
};

#endif // CELL_H
