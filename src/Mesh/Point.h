#ifndef POINT_H
#define POINT_H

#include "numvector.h"
#include <memory>

///
/// 2D Point
///

class Point : public numvector<double, 2>
{

public:

    /// Default constructor
    Point(double val = 0.0) : numvector(val) {}

    /// Construct with defined values
    Point(const numvector<double, 2>& coord) : numvector(coord) {}

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

    /// check if point is equal to given
    bool isEqual(const Point& p) const { return true ? (r[0] == p.x() && r[1] == p.y()) : false; }

    /// check if point is equal to given
    bool isEqual(const std::shared_ptr<Point>& p) const { return isEqual(*p); }

    /// save global node number just for VTK export simplification
    int gNumber;
};

#endif // POINT_H
