#ifndef PATCH_H
#define PATCH_H

#include <string>
#include <vector>
#include "EdgeBoundary.h"

///
/// Group og edges on the flow domain boundary where some boundary condition should be applied
///

class Patch
{

public:

    /// Name of patch
    std::string patchName;

    /// List of edges on the patch
    std::vector<std::shared_ptr<EdgeBoundary>> edgeGroup;

public:

    /// Default constructor
    Patch() {}

    /// Destructor
    ~Patch() {}
};

#endif // PATCH_H
