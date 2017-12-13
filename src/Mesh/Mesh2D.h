#ifndef MESH2D_H
#define MESH2D_H

#include "numvector.h"
#include "Point.h"
#include "Cell.h"
#include "FluxLLF.h"
#include "EdgeInternal.h"
#include "EdgeBoundaryInfty.h"

#include <fstream>


class Mesh2D
{

private:

    //- File ofstream for mesh export
    std::ofstream writer;

public:

    //- Space step in x direction
    //double hx;

    //- Space step in y direction
    //double hy;

    //- Number of internal cells  //?
    int nCells;


    //- Coordinates of nodes (x,y)
    std::vector<Point> nodes;

    //- Horizontal edges (nodeLeft, nodeRight)
    std::vector<Edge*> edgesHor;

    //- Vertical edges (nodeDown, nodeUp)
    std::vector<Edge*> edgesVer;

    //- Boundary edges
    std::vector<EdgeBoundaryInfty*> edgesBound;

    //- Mesh cells (edgeDown, edgeUp, edgeLeft, edgeRight)
    std::vector<Cell*> cells;

public:

    //- Default constructor
    Mesh2D() {}

    //- Construct uniform rectangular mesh by number of cells and size of flow domain
    Mesh2D(int nx, int ny, double Lx, double Ly);

    //- Destructor
    ~Mesh2D();

    //- Export mesh
    void exportMesh();

};// end Mesh 2D

#endif // MESH2D_H

