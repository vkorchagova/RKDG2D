

#include "numvector.h"
#include "defs.h"
#include "Params.h"
#include "Cell.h"


///
/// Rotation
///

//- rotate coordinate system clockwise
numvector<double, dimPh> rotate(const numvector<double, dimPh>& sol, const Point& n);

//- rotate coordinate system counter-clockwise
numvector<double, dimPh> inverseRotate(const numvector<double, dimPh>& sol, const Point& n);


///
/// Integration
///

//- 2D Gauss integration of scalar function
double integrate( const Cell& cell, const std::function<double(const Point &)>& f);

//- 2D Gauss integration of vector function
numvector<double,5> integrate( const Cell& cell, const std::function<numvector<double, 5>(const Point&)>& f);


