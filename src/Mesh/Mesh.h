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

/// Proc rank 
extern int myRank;

/// Number of procs
extern int numProcsTotal;

/// Status
extern MPI_Status status;

/// Debug
extern bool debug;

/// Log file to save data
extern std::ofstream logger;

///
/// Mesh
///

class Mesh
{

private:

    /// Reference to MPI buffers
    Buffers& buf;

    /// Find neighbours for given cell
    void findNeighbourCells(const std::shared_ptr<Cell>& cell);

    /// Find vertex-neighbours for given cell
    void findNeighbourCellsVertex(const std::shared_ptr<Cell>& cell);

    /// Add cell to proc patch
    void addToProcPatch(const std::shared_ptr<Cell>& cell, int numProc);

    void createPhysicalPatch(const std::vector<std::shared_ptr<Edge>>& edgeGroup, const std::string& pName);

    /// renumeration of edges: bound|inner|proc
    void structurizeEdges();

    /// running for neighbor cells by vertices
    bool moveCell(const std::shared_ptr<Cell>& cell, std::shared_ptr<Cell>& startCell, std::vector<std::shared_ptr<Cell>>& neibCellVertexOneEdge);
   // bool moveCell(const std::shared_ptr<Cell>& cell, std::shared_ptr<Cell>& startCell, std::shared_ptr<Cell>& nextCell, bool& ccw);

    /// filter neighbor cells to get convex stencil
    void filterNeibCellsVertex(const std::shared_ptr<Cell>& cell);

public:	

    /// Nodes of the mesh
    std::vector<std::shared_ptr<Point>> nodes;

    /// Edges
    std::vector<std::shared_ptr<Edge>> edges;

    /// Mesh cells (edge1, ..., edgek) CCW
    std::vector<std::shared_ptr<Cell>> cells;

    /// Groups of physical boundary edges
    std::vector<Patch> patches;

    /// Groups of shadow cells between neib procs
    std::vector<ProcPatch> procPatches;

    /// Global numeration of nodes
    std::vector<int> globalNodeNumber;

    /// Global numeration of edges
    std::vector<int> globalEdgeNumber;

    /// Global numeration of cells
    std::vector<int> globalCellNumber;

    /// Group of cells special for finDiff additional limitation
    std::vector<int> finDiffGroup;

    /// Number of nodes
	int nNodes;

    /// Numbers of edges
    int nRealEdges;

    /// MAY BE DEPRECATED Total number of edges (includes edges for neib proc cells)
    int nTotalEdges;

    /// Numbers of real cells
    int nRealCells;

    /// Number of ghost cells
    int nEdgesBound;

    /// Number of patches
    int nPatches; //needs to be delete? see applyBoundary function

    /// Number of neighbor processors
    int nNeibProcs;

    /// Number of cells in full mesh
    int nCellsGlob;

public:

    /// Constructor
    Mesh(Buffers& buf) : buf(buf) {};

    /// Constructor
    Mesh(std::string& fileName, Buffers& buf);

    /// Destructor
    ~Mesh();

    /// Import mesh from file "mesh2D"
    void importMesh(std::string& fileName);

    /// Export VTK for vertices too
    void exportMeshVTK_polyvertices(std::ostream& writer) const;

};// end Mesh

#endif // MESH_H

