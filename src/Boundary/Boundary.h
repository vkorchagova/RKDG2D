#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "numvector.h"
#include "Point.h"
#include <string>


class Boundary
{
public:

    //- Type of BC
    std::string type;

    //- Default constructor
    Boundary() : type("not implemented") {};

    //- Apply boundary conditions
    virtual numvector<double, 5> applyBoundary(const numvector<double, 5>& sol = {0.0, 0.0, 0.0, 0.0, 0.0}, const Point& n = Point({0.0,0.0})) const = 0;
};

#endif // BOUNDARY_H
