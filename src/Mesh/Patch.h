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
    std::string patchName;

    //- List of numbers of cells in patch
    std::vector<int> cellGroup;

    //- Pointer to boundary condition
    //Boundary bc;

    //- Default constructor
    Patch() {}

    //- Destructor
    ~Patch() {}
};

#endif // PATCH_H
