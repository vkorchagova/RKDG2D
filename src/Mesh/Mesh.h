#ifndef MESH_H
#define MESH_H

#include "numvector.h"
#include "Params.h"
#include "defs.h"
#include "compService.h"
#include "Edge.h"
#include "Patch.h"
#include "Cell.h"

#include <vector>
#include <fstream>
#include <memory>


class Mesh
{

private:

    //- Find neighbours for given cell
    void findNeighbourCells (const std::shared_ptr<Cell>& cell) const;

public:	

    //- Nodes of the mesh
    std::vector<Point> nodes;

    //- Edges
    std::vector<Edge> edges;

    //- Mesh cells (edge1, ..., edgek) CCW
    std::vector<Cell> cells;

    //- Groups of boundary edges
    std::vector<Patch> patches;

    //- Global numeration of nodes
    std::vector<int> globalNodeNumber;

    //- Global numeration of edges
    std::vector<int> globalEdgeNumber;

    //- Global numeration of cells
    std::vector<int> globalCellNumber;

    //- Number of nodes
	int nNodes;

    //- Numbers of edges
    int nEdges;

    //- Numbers of real cells
    int nRealCells;

    //- Nuber of ghost cells
    int nGhostCells;

    //- Number of patches
    int nPatches; //needs to be delete? see applyBoundary function

public:

    //- Constructor
    Mesh(std::string& fileName);

    //- Destructor
    ~Mesh();

    //- Import mesh from file "mesh2D"
    void importMesh(std::string& fileName);

    //- Export VTK only for cells
    void exportMeshVTK(std::ostream& writer) const;

    //- Export VTK for vertices too
    void exportMeshVTK_polyvertices(std::ostream& writer) const;

};// end Mesh

#endif // MESH_H

