#ifndef GAUSSINTEGRATOR_H_
#define GAUSSINTEGRATOR_H_

#include "numvector.h"
#include <functional>
#include <math.h>

namespace std
{

class GaussIntegrator
{
private:

    //- gauss points
    numvector<double,2> gPoints1D;
    numvector<numvector<double,2>,4> gPoints2D;

    //- weights
    numvector<double,2> gWeights1D;
    numvector<double,4> gWeights2D;

    //- local [-1;1] to global [a;b]
    double localToGlobal(double x, double a, double b);

    //- local [-1,1]x[-1,1] to global rectangular cell
    numvector<double, 2> localToGlobal(numvector<double, 2> coord, numvector<numvector<double,2>,4> nodes);

    //- compute area of rectangle
    double getArea(numvector<numvector<double,2>,4> nodes)
        {return (nodes[1][0] - nodes[0][0])*(nodes[2][1] - nodes[1][1]);}

public:
    GaussIntegrator();
    ~GaussIntegrator(){}

    //- 1D Gauss integration of scalar function
    //double integrate( const function<double(double)>& f, double a, double b);

    //- 1D Gauss integration of vector function
    numvector<double,5> integrate( const function<numvector<double,5>(double)>& f, double a, double b);

    //- 2D Gauss integration (only for rectangles defined by nodes)
    //double integrate( const function<double(const numvector<double, 2>&)>& f, const numvector<numvector<double,2>,4>& nodes);

    //- 2D Gauss integration of vector function (only for rectangles defined by nodes)
    numvector<double,5> integrate( const function<numvector<double,5>(const numvector<double, 2>&)>& f, const numvector<numvector<double,2>,4>& nodes);

};// End GaussIntegrator

} // End namespace std

#endif // GAUSSINTEGRATOR_H
