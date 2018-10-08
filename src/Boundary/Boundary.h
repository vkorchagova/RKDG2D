#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "numvector.h"
#include "Point.h"
#include <string>

///
/// Abstract class for boundary conditions
///

class Boundary
{
public:

    /// Type of BC
    ///
    /// Default value: "not implemented"
    std::string type;

    /// Constructor
    Boundary() : type("not implemented") {};

    /// Apply boundary conditions
    ///
    /// This function is used for computation of numerical fluxes to get the value of solution outside boundary cell
    ///
    /// @param  sol solution near the point n inside boundary cell
    /// @param  n   Normal direction related to boundary edge
    ///
    /// @return the value of solution outside boundary cell    
    virtual numvector<double, 5> applyBoundary(const numvector<double, 5>& sol = {0.0, 0.0, 0.0, 0.0, 0.0}, const Point& n = Point({0.0,0.0})) const = 0;
};

#endif // BOUNDARY_H
