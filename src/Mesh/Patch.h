/// -----------------------------------
/// Patch class of egde boundary groups
/// -----------------------------------

#ifndef PATCH_H
#define PATCH_H

#include <string>
#include <vector>
#include "Edge.h"

// Default patch
class Patch
{

public:

    //- Name of patch
    std::string name;

    std::vector<std::shared_ptr<Edge>> edgeGroup;

    //- List of numbers of cells in patch
    std::vector<std::shared_ptr<Cell>> cellGroup;

    //- Default constructor
    Patch() {}

    //- Construct by name
    Patch(std::string pN) : name(pN) {}

    //- Destructor
    ~Patch() {}
};


// MPI proc patch
class ProcPatch : public Patch
{

public:

    //- Proc num
    int procNum;

    //- Constructor
    ProcPatch(std::string pN, int pNum) : Patch(pN), procNum(pNum) {}

    //- Destructor
    ~ProcPatch() {}

};

#endif // PATCH_H
