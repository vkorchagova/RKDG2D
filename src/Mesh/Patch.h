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

    //- List of numbers of cells in patch
    std::vector<std::shared_ptr<Cell>> cellGroup;

    //- Pointer to boundary condition
    //Boundary bc;

    //- Default constructor
    Patch() {}

    //- Construct by name
    Patch(std::string pN) : name(pN) {}

    //- Destructor
    ~Patch() {}
};

#endif // PATCH_H
