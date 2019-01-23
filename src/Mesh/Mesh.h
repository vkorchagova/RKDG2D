#ifndef MESH_H
#define MESH_H

#include "numvector.h"
#include "Params.h"
#include "Edge.h"
#include "Patch.h"
#include "Cell.h"

#include "Buffers.h"

#include <vector>
#include <fstream>
#include <memory>

//- Proc rank 
extern int myRank;

//- Size of ...
extern int numProcsTotal;

//- Status
extern MPI_Status status;

//- Request
extern MPI_Request request;


class Mesh
{

private:

    //- Reference to MPI buffers
    Buffers& buf;

    //- Find neighbours for given cell
    void findNeighbourCells(const std::shared_ptr<Cell>& cell);

    //- Add cell to proc patch
    void addToProcPatch(const std::shared_ptr<Cell>& cell, int numProc);

    void createPhysicalPatch(const std::vector<std::shared_ptr<Edge>>& edgeGroup, const std::string& pName);

    std::shared_ptr<Cell> makeGhostCell(const std::shared_ptr<Edge>& e);

    void structurizeEdges();

public:	

    //- Nodes of the mesh
    std::vector<std::shared_ptr<Point>> nodes;

    //- Edges
    std::vector<std::shared_ptr<Edge>> edges;

    //- Mesh cells (edge1, ..., edgek) CCW
    std::vector<std::shared_ptr<Cell>> cells;

    //- Groups of physical boundary edges
    std::vector<Patch> patches;

    //- Groups of shadow cells between neib procs
    std::vector<ProcPatch> procPatches;

    //- Global numeration of nodes
    std::vector<int> globalNodeNumber;

    //- Global numeration of edges
    std::vector<int> globalEdgeNumber;

    //- Global numeration of cells
    std::vector<int> globalCellNumber;

    //- Number of nodes
	int nNodes;

    //- Numbers of edges
    int nRealEdges;

    //- Total number of edges (includes edges for neib proc cells)
    int nTotalEdges;

    //- Numbers of real cells
    int nRealCells;

    //- Number of ghost cells
    int nGhostCells;

    //- Number of patches
    int nPatches; //needs to be delete? see applyBoundary function

    //- Number of neighbor processors
    int nNeibProcs;

    //- Number of cells in full mesh
    int nCellsGlob;

    // //- Full global numeration of cells
    // std::vector<int> fullGlobalMap;

    // //- Array of displacements for MPI_Gatherv
    // std::vector<int> mpiDispl;

    // //- nRealCells on each process
    // std::vector<int> nCellsPerProc;

public:

    //- Constructor
    Mesh(std::string& fileName, Buffers& buf);

    //- Destructor
    ~Mesh();

    //- Import mesh from file "mesh2D"
    void importMesh(std::string& fileName);

    //- Export VTK for vertices too
    void exportMeshVTK_polyvertices(std::ostream& writer) const;

};// end Mesh

#endif // MESH_H

