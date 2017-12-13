#ifndef CELL_H
#define CELL_H

/// ------------------------------
/// Cell class
/// ------------------------------
/// Consists of 4 edges
/// Know its Gauss points
/// Can evaluate RHS inside cell
/// ------------------------------

#include "Point.h"
//#include "Edge.h"
#include "numvector.h"
#include "Problem.h"

#include <functional>
#include <math.h>

class Edge;
//class Problem;

class Cell
{

private:

    //- Number of Gauss points
    static const int nGP = 4;

    //- Gauss points
    numvector<Point, nGP> gPoints2D;

    //- Gauss weights
    numvector<double, nGP> gWeights2D;

private:
    //- Area of cell
    double area;

    //- Mass center of cell
    Point center;

    Point hxhy;

    //- Check if point belongs cell
    bool insideCell(const Point& point) const;

    //- local [-1,1]x[-1,1] to global rectangular cell
    Point localToGlobal(const Point& point) const;

    //- Set Gauss points
    void setGaussPoints();

    //- Set basis functions
    void setBasisFunctions();
    
    //- Default constructor
    Cell() {};


public:

    /// geometric variables

    //- Compute hx, hy
    const Point& h() const { return hxhy; };

    //- Number of cell
    //debug
    int number;

    //- Number of edges
    static const int nEdges = 4;

    //- Edges define cell
    numvector<Edge*, nEdges> edges;

    //- Compute area of cell
    double getArea() const { return area; };

    //- Calculate center of cell
    const Point& getCellCenter() const { return center; };

    /// RKDG variables

    //- Pointer to problem
    const Problem* problem;

    //- List of basis functions
    //std::function<double(const Point&)> phi[nShapes];
    std::vector<std::function<double(const Point&)>> phi;

    //- Gradient of basis functions
    std::vector<std::function<Point(const Point&)>> gradPhi;

public:


    //- Construct cell using numvector of edges
    Cell(const numvector<Edge*, 4>& edges);

    //- Destructor
    ~Cell() {};

    /// geometric methods

    //- Calculate coordinates of cell nodes
    numvector<const Point*, nEdges> getCellCoordinates() const;


    /// RKDG methods

    //- Set problem
    void setProblem(const Problem& prb);

    //- Reconstruct solution
    numvector<double, 5> reconstructSolution(const Point& point) const;

    //- Get coefficients for initial conditions
    numvector<double, 5 * nShapes> getLocalInitialConditions(std::function<numvector<double,5>(const Point& point)>& init) const;

    //- Calculate local RHS
    //numvector<double, 5 * nShapes> getLocalRHS() const;

    //- 2D Gauss integration of vector function
    numvector<double,5> integrate( const std::function<numvector<double, 5>(const Point&)>& f) const;


};

#endif // CELL_H
