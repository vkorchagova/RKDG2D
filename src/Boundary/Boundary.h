#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "numvector.h"


class Boundary
{
public:

    //- Default constructor
    Boundary() {};

    //- Apply boundary conditions
    virtual numvector<double, 5> applyBoundary(const numvector<double, 5>& sol = {0.0, 0.0, 0.0, 0.0, 0.0}) const = 0;
};

#endif // BOUNDARY_H
