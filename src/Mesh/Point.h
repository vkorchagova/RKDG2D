/// ------------------------------
/// 2D Point class
/// ------------------------------


#ifndef POINT_H
#define POINT_H

#include "numvector.h"

class Point
{

private:

    numvector<double, 2> coord;

public:

    //- Default constructor
    Point() { coord[0] = 0.0; coord[1] = 0.0; }

    //- Construct with defined values
    Point(double x, double y) { set(x,y); }

    //- Destructor
    ~Point() {}

    //- Get x coordinate
    double x() const { return coord[0]; }

    //- Get y coordinate
    double y() const { return coord[1]; }

    //- Set x, y values
    void set(double x, double y) { coord[0] = x; coord[1] = y; }

    //- Overloaded "+=" operator
    Point& operator+=(Point& p) { coord[0] += p.coord[0]; coord[1] += p.coord[1]; return *this; }
};

#endif // POINT_H
