#ifndef PATCH_H
#define PATCH_H

#include <string>
#include <vector>
#include "Edge.h"

///
/// Geometrical part of boundary where one boundary condition should be apply
///

class Patch
{

public:

    /// Name of patch
    std::string name;

    /// List of edges the patch consists of
    std::vector<std::shared_ptr<Edge>> edgeGroup;

    /// List of cells adjacent to patch
    std::vector<std::shared_ptr<Cell>> cellGroup;

    /// Default constructor
    Patch() {}

    /// Construct by name
    Patch(std::string pN) : name(pN) {}

    /// Destructor
    ~Patch() {}
};


///
/// Special part of boundary which placed between two parts of mesh processed by different procs
///

class ProcPatch : public Patch
{

public:

    /// Proc num
    int procNum;

    /// Inner cell group
    std::vector<std::shared_ptr<Cell>> innerCellGroup;

    /// Constructor
    ProcPatch(std::string pN, int pNum) : Patch(pN), procNum(pNum) {}

    /// Destructor
    ~ProcPatch() {}

};

#endif // PATCH_H
