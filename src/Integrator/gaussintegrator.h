#ifndef GAUSSINTEGRATOR_H_
#define GAUSSINTEGRATOR_H_

#include "numvector.h"
#include <functional>
#include <math.h>
#include "Cell.h"
#include "Point.h"

class Cell;


class GaussIntegrator
{
private:

    //- gauss points
    Point gPoints1D;
    numvector<Point,4> gPoints2D;

    //- weights
    numvector<double,2> gWeights1D;
    numvector<double,4> gWeights2D;

    //- local [-1;1] to global [a;b]
    double localToGlobal(double x, double a, double b);

    //- local [-1,1]x[-1,1] to global rectangular cell
    Point& localToGlobal(Point& point, Cell& cell);

public:
    GaussIntegrator();
    ~GaussIntegrator(){}

    //- 1D Gauss integration of scalar function
    //double integrate( const function<double(double)>& f, double a, double b);

    //- 1D Gauss integration of vector function
    numvector<double,5> integrate( const std::function<numvector<double,5>(double)>& f, double a, double b);

    //- 2D Gauss integration (only for rectangles defined by nodes)
    //double integrate( const function<double(const numvector<double, 2>&)>& f, const numvector<numvector<double,2>,4>& nodes);

    //- 2D Gauss integration of vector function
    numvector<double,5> integrate( const std::function<numvector<double,5>(const Point&)>& f, Cell& cell);

};// End GaussIntegrator


#endif // GAUSSINTEGRATOR_H
