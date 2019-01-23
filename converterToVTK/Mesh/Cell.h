#ifndef CELL_H
#define CELL_H

/// ------------------------------
/// Cell class
/// ------------------------------
/// Consists of edges
/// Know its Gauss points
/// ------------------------------


#include "numvector.h"

#include <functional>
#include <memory>
#include <math.h>
#include <algorithm>

#include "Point.h"
#include "Edge.h"

class Edge;

class Cell
{

public:

    /// geometric variables

    int number; 
    
    //- Number of Gauss points
    int nGP;

    //- Gauss points
    std::vector<Point> gPoints2D;
    
    //- Gauss weights
    std::vector<double> gWeights2D;
    
    //- Jacobian
    std::vector<double> J;
    
    //- Area of cell
    double area;
    
    //- Mass center of cell
    Point center;

    //- Gramian matrix
    std::vector<std::vector<double>> gramian;

    //- local [-1,1]x[-1,1] to global rectangular cell
    Point localToGlobal(const Point& localPoint) const;

    //- Nodes in cell
    //std::vector<std::shared_ptr<Point>> nodes;
    std::vector<std::shared_ptr<Point>> nodes;

    //- Edges in cell
    //std::vector<std::shared_ptr<Edge>> edges;
    std::vector<std::shared_ptr<Edge>> edges;

    //- Number of entities (nodes or edges)
    int nEntities;

    //- Neighbour cells
    //std::vector<std::shared_ptr<Cell>> neibCells;
    std::vector<std::shared_ptr<Cell>> neibCells;

    //- Return area of cell
    double getArea() const { return area; }

    //- Return center of cell
    const Point& getCellCenter() const { return center; }

    //- Calculate element area
    void setArea();

    //- Set cell center
    void setCellCenter();

    //- Define Jacobian function
    void setJacobian();

    //- Set Gauss points
    void setGaussPoints();


    //- Construct cell using vectors of nodes and edges
    //Cell(const std::vector<std::shared_ptr<Point>> &nodes, const std::vector<std::shared_ptr<Edge>> &edges);
    Cell(const std::vector<std::shared_ptr<Point>> &nodes, const std::vector<std::shared_ptr<Edge>> &edges);
    
    //- Copy constructor
    //Cell(const Cell&) = delete;

    //- Destructor
    ~Cell() {}

    /// geometric methods

    //- Calculate coordinates of cell nodes
    std::vector<Point> getCellCoordinates() const;

    //- Find neighbour cells
    void findNeighbourCells();
};

#endif // CELL_H
