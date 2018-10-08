#ifndef POINT_H
#define POINT_H

#include "numvector.h"

///
/// Geometric point
///

class Point : public numvector<double, 2>
{

public:

    /// Default constructor
    Point(double val = 0.0) : numvector(val) {}

    /// Construct with defined values
    Point(const numvector<double, 2>& coord) : numvector(coord) {}

    /// Copy constructor
    Point(const Point& p) = default;

    /// Destructor
    ~Point() {}

    /// Get x coordinate
    const double& x() const { return r[0]; }

    /// Set x coordinate
    double& x() { return r[0]; }

    /// Get y coordinate
    const double& y() const { return r[1]; } 

    /// Set y coordinate
    double& y() { return r[1]; }
};

///
/// Mesh node
///

class Node : public Point
{

public:

    /// Number of node
    int number;

    /// Constructor
    Node(const Point& p, int num) : Point(p), number(num ) {}

};


#endif // POINT_H
