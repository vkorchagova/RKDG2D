#ifndef CELL_H
#define CELL_H

/// ------------------------------
/// Cell class
/// ------------------------------
/// Consists of 4 edges
/// Know its Gauss points
/// Can evaluate RHS inside cell
/// ------------------------------

#include "Edge.h"
#include "numvector.h"
#include "Problem.h"
#include "Point.h"
#include <functional>

//class Edge;
class Problem;


class Cell
{

private:

    //- Gauss points
    numvector<Point,4> gPoints2D;

    //- Gauss weights
    numvector<double,4> gWeights2D;

    //- Number of Gauss points
    int nGP;



private:

    //- Number of basis functions
    static const int nShapes = 3;

    //- Compute hx, hy
    void getSteps();

    //- Compute area of cell
    void getArea();

    //- Calculate center of cell
    void getCellCenter();

    //- Check if point belongs cell
    bool insideCell(Point& point);

    //- local [-1,1]x[-1,1] to global rectangular cell
    Point localToGlobal(Point& point);

    //- Set Gauss points
    void setGaussPoints();

    //- Set basis functions
    void setBasisFunctions();


public:

    /// geometric variables


    //- Space step in x direction
    double hx;

    //- Space step in y direction
    double hy;

    //- Number of cell
    int number;

    //- Number of edges
    static const int nEdges = 4;

    //- Edges define cell
    numvector<Edge*, nEdges> edges;

    //- Area of cell
    double area;

    //- Mass center of cell
    Point center;


    /// RKDG variables

    //- Pointer to problem
    Problem* problem;

    //- Number of 2D Gauss points for cell
    //static const int nGP = 4;

    //- List of basis functions
    //std::function<double(const Point&)> phi[nShapes];
    std::vector<std::function<double(const Point&)>> phi;

    //- Gradient of basis functions
    std::function<Point(const Point&)> gradPhi[nShapes];

    //- Solution coeffs on cells on previous time step
    numvector<double, 5 * nShapes> alphaPrev;

    //- Solution coeffs on cells on next time step
    numvector<double, 5 * nShapes> alphaNext;

public:

    //- Default constructor
    Cell();

    //- Construct cell using numvector of edges
    Cell(numvector<Edge *,4> &edges);

    //- Destructor
    ~Cell();

    //- Copy constructor
    Cell(const Cell& c);

    //- Overload of "=" operator
    Cell& operator=(const Cell& c);


    /// geometric methods

    //- Calculate coordinates of cell nodes
    numvector<Point*, nEdges> getCellCoordinates();


    /// RKDG methods

    //- Set problem
    void setProblem(Problem& prb);

    //- Reconstruct solution
    numvector<double, 5> reconstructSolution(const numvector<double, nShapes * 5>& alpha, Point& point);

    //- Set initial conditions
    void setLocalInitialConditions(std::function<numvector<double,5>(const Point& point)>& init);

    //- Calculate local RHS
    numvector<double, 5 * nShapes> getLocalRHS();

    //- 2D Gauss integration of vector function
    numvector<double,5> integrate( const std::function<numvector<double,5>( Point&)>& f);


};

#endif // CELL_H
