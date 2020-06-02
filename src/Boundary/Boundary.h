#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "numvector.h"
#include "Patch.h"
#include "Params.h"
#include <string>

/// 
/// Abstract class for boundary condition
/// 

class Boundary
{
public:

    /// Type of boundary condition
    std::string type;

    /// Constant reference to geometrical boundary
    const Patch& patch;

    /// Default constructor
    Boundary(const Patch& p) : type("not implemented"), patch(p) {};

    /// Get solution in outer side (pure virtual)
    virtual numvector<double, dimPh> getSolOuter(
        const numvector<double, dimPh>& solInner, 
        const Point& n = Point({0.0,0.0})
    ) const = 0;

};

#endif // BOUNDARY_H
