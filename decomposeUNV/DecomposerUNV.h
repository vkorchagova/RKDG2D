#ifndef DECOMPOSERUNV_H
#define DECOMPOSERUNV_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>

/// Mesh Ideas UNV decomposer based on METIS program
///
/// Result of partition writes to sequence of files looks like
///     mesh2D.0
///     mesh2D.1
///     ...
///
/// File format for output is the special format for RKDG code
/// - see more after class definition
///
/// IDEAS UNV mesh: the sequence of blocks like:
///     -1
/// universal dataset number according to IDEAS notation
/// (see more: http://www.sdrl.uc.edu/sdrl/referenceinfo/universalfileformats/file-format-storehouse)
/// content of section
///     -1
///
/// More about METIS:
///     http://glaros.dtc.umn.edu/gkhome/metis/metis/overview
///


//- Boundary between UNV sections
const std::string SEPARATOR("    -1");


class DecomposerUNV
{

private:

    /// Variables

    //- fstream for UNV mesh
    std::ifstream reader;

    //- fstream for RKDG mesh
    std::ofstream writer;

    //- nodes
    std::vector<std::vector<double>> nodes;

    //- edges = numbers of nodes
    std::vector<std::vector<int>> edges;

    //- cells = numbers of edges
    std::vector<std::vector<int>> cellsAsEdges;
    std::vector<std::vector<int>> cellsAsNodes;

    //- cell centers
    std::vector<std::vector<double>> cellCenters;

    //- adjoint cells for edges
    std::vector<std::vector<int>> adjEdgeCells;

    //- edge normals
    std::vector<std::vector<double>> edgeNormals;

    //- patch names
    std::vector<std::string> patchNames;

    //- patch edge groups
    std::vector<std::vector<int>> patchEdgeGroups;
    
    //- partition to subregions
    std::vector<int> partition;


    /// Methods

    //- Read number of section
    int readTag();

    //- Skip section
    void skipSection();

    //- Read nodes
    void readNodes();

    //- Read elements (edges + cells)
    void readElements();

    //- Read patches
    void readPatches();

    //- Find edges for cell defined by nodes
    void getElementEdges(const std::vector<int>& nodeNumbers);

    //- Find center of cell vertices
    void getCellCenter(const std::vector<int>& nodes);

    //- Find normal vectors for edges
    void getEgdeNormals();

    //- Find ajoint cells for each edge
    void setAdjointCells();


public:

    //- Constructor
    DecomposerUNV(std::string unvMeshFile, std::string rkdgMeshFile);

    //- Destructor
    ~DecomposerUNV();

    //- Import UNV mesh
    void importUNV ();

    //- Export RKDG mesh
    void exportRKDG();
    
    //- Export mesh elements for METIS utility
    void exportMETIS() const;
    
    //- Read results of partition
    void importPartition(std::string metisPartFile);
    
    //- Export VTK visualisation
    void exportVTK() const;
    
    //- Export part of mesh
    void exportPartMeshRKDG(int numDom) const;
    
};


#endif // DECOMPOSERUNV_H

/// RKDG mesh format for decomposed mesh
/// ------------------------------------
///
/// Each number are related to appropriate part of mesh
///
/// $Nodes
///     number_of_nodes
///     global_number x y
///     ...
/// $EndNodes
/// $Edges //at first - boundary edges, after - internal
///     number_of_boundary_edges
///     total_number_of_edges
///     global_number is_edge_on_boundary node_1 node_2
///     ...
/// $EndEdges
/// $Cells
///     number_of_cells
///     global_number number_of_edges_in_cell node1 node2 ... edge1 edge2 ... // number_of_nodes = number_of_edges
///     ...
/// $EndCells
/// $AdjointCellsForEdges
///     total_number_of_edges
///     cell1 cell2
///     ...
/// $EndAdjointCellsForEdges
/// $EdgeNormals
///     total_number_of_edges
///     n_x n_y
/// /// ...
/// $EndEdgeNormals
/// $Patches
///     number_of_patches
///     patch_name1
///         number_of_edges_in_patch
///         edgenum1
///         edgenum2
///         ...
///     patch_name_2
///         ...
/// $EndPatches
