#pragma once

#include "numvector.h"
#include "Edge.h"
#include "Cell.h"
#include <fstream>

class Point;
class Cell;
class Edge;


class Mesh2D
{
private:

    //- File ofstream for mesh export
    std::ofstream writer;

public:

	//- Space step in x direction
	double hx;

	//- Space step in y direction
	double hy;

    //- Number of internal cells
    int nInternalCells;

    //- Number of ghost cells
    //int nGhostCells;



    //- Coordinates of nodes (x,y)
    std::vector<Point> nodes;

    //- Horizontal edges (nodeLeft, nodeRight)
    std::vector<Edge> edgesHor;

    //- Vertical edges (nodeDown, nodeUp)
    std::vector<Edge> edgesVer;

    //- Mesh cells (edgeDown, edgeUp, edgeLeft, edgeRight)
    std::vector<Cell> cells;

    //- Internal mesh cells
    std::vector<Cell> internalCells;

    //- Ghost mesh cells (for boundary)
    //std::vector<Cell*> ghostCells;

public:

    //- Default constructor
    Mesh2D() {}

    //- Construct uniform rectangular mesh by number of cells and size of flow domain
	Mesh2D(int nx, int ny, double Lx, double Ly);

    //- Destructor
    ~Mesh2D();


    //- Get coordinates of nodes of given mesh cell
	numvector<numvector<double,2>,4> getCellCoordinates(int iCell);

    //- Export mesh
    void exportMesh();

};// end Mesh 2D


