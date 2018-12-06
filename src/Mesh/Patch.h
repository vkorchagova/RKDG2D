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

    //- List of edges in the patch
    std::vector<std::shared_ptr<Edge>> edgeGroup;

    //- Default constructor
    Patch() {}

    //- Destructor
    ~Patch() {}
};

#endif // PATCH_H
