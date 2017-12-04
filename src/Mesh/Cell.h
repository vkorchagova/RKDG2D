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

class Edge;
class Problem;

class Cell
{

private:

    //- Space step in x direction
    double hx;

    //- Space step in y direction
    double hy;

private:

    //- Compute hx, hy
    void getSteps();

    //- Compute area of cell
    void getArea();

    //- Calculate center of cell
    void getCellCenter();

public:

    /// geometric variables

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
    static const int nGP = 4;

public:

    //- Default constructor
    Cell();

    //- Construct cell using numvector of edges
    Cell(numvector<Edge *,4> &edges);

    //- Destructor
    ~Cell();


    /// geometric methods

    //- Calculate coordinates of cell nodes
    numvector<Point*, nEdges> getCellCoordinates();




    /// RKDG methods

    //- Calculate local RHS
    double getLocalRHS();


};

#endif // CELL_H
