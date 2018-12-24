/// -----------------------------------
/// Patch class of egde boundary groups
/// -----------------------------------

#ifndef PATCH_H
#define PATCH_H

#include <string>
#include <vector>
#include "Edge.h"

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

    //- Apply boundary
    // void applyBoundary();
};

#endif // PATCH_H
