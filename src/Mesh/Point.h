/// ------------------------------
/// 2D Point class
/// ------------------------------

#ifndef POINT_H
#define POINT_H

#include "numvector.h"


class Point : public numvector<double, 2>
{

public:

    //- Number of node
    int number;

    //- Default constructor
    Point(double val = 0.0) : numvector(val) {}

    //- Construct with defined values
    Point(const numvector<double, 2>& coord) : numvector(coord) {}

    //- Copy constructor
    Point(const Point& p) = default;

    //- Destructor
    ~Point() {}

    //- Get/Set x coordinate
    const double& x() const { return r[0]; }
    double& x() { return r[0]; }

    //- Get/Set y coordinate
    const double& y() const { return r[1]; }
    double& y() { return r[1]; }
};

#endif // POINT_H
