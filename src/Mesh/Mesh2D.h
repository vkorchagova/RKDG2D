#ifndef MESH2D_H
#define MESH2D_H

#include "numvector.h"
#include "Point.h"
#include "Patch.h"
#include "Cell.h"
#include "FluxLLF.h"
#include "EdgeInternal.h"
#include "EdgeBoundaryDiagProjection.h"

#include <fstream>
#include <memory>

///
/// Two-dimensional mesh
///

class Mesh2D
{

private:

    /// File ofstream for mesh export
    mutable std::ofstream writer;

    /// File ifstream for mesh import
    mutable std::ifstream reader;

    /// Find neighbours for given cell
    void findNeighbourCells (const std::shared_ptr<Cell>& cell) const;

public:

    /// Number of cells
    int nCells;

    /// Number of edges
    int nEdges;

    /// Number of edges on boundary of flow domain
    int nBoundEdges;

    /// Sum of entities in each cell (need for VTK writer)
    int nEntitiesTotal;

    /// Coordinates of nodes (x,y)
    std::vector<Node> nodes;

    /// Edges
    std::vector<std::shared_ptr<Edge>> edges;

    /// Mesh cells (edge_1, ..., edge_k counterclockwise)
    std::vector<std::shared_ptr<Cell>> cells;

    /// Groups of boundary edges
    std::vector<Patch> patches;

public:

    /// Construct mesh by import from UNV file
    Mesh2D(std::string fileName, const Problem& prb);

    /// Destructor
    ~Mesh2D();

    /// Import mesh
    void importMesh(std::string fileName, const Problem &prb);

    /// Export arbitrary 2D mesh in custom RKDG format like .msh
    void exportMesh() const;

    /// Export VTK only with cell centers
    void exportMeshVTK(std::ostream& writer) const;

    /// Export VTK with nodes and cell centers
    void exportMeshVTK_polyvertices(std::ostream& writer) const;

};// end Mesh 2D

#endif // MESH2D_H

