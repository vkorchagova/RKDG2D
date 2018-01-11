#ifndef MESH2D_H
#define MESH2D_H

#include "numvector.h"
#include "Point.h"
#include "Cell.h"
#include "FluxLLF.h"
#include "EdgeInternal.h"
#include "EdgeBoundaryInfty.h"
#include "EdgeBoundarySlip.h"
#include "EdgeBoundaryOpen.h"

#include <fstream>
#include <memory>


class Mesh2D
{

private:

    //- File ofstream for mesh export
    mutable std::ofstream writer;

public:

    //- Number of internal cells  //?
    int nCells;

    //- Coordinates of nodes (x,y)
    std::vector<Point> nodes;

    //- Horizontal edges (nodeLeft, nodeRight)
    std::vector<std::shared_ptr<Edge>> edgesHor;

    //- Vertical edges (nodeDown, nodeUp)
    std::vector<std::shared_ptr<Edge>> edgesVer;

    //- Mesh cells (edgeDown, edgeUp, edgeLeft, edgeRight)
    std::vector<std::shared_ptr<Cell>> cells;

public:

    //- Default constructor
    // Mesh2D() {}

    //- Construct uniform rectangular mesh by number of cells and size of flow domain
    Mesh2D(int nx, int ny, double Lx, double Ly);

    //- Destructor
    ~Mesh2D();

    //- Export mesh
    void exportMesh() const;

};// end Mesh 2D

#endif // MESH2D_H

