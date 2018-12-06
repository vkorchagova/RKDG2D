#ifndef COMPSERVICE_H_
#define COMPSERVICE_H_

#include "numvector.h"
#include "defs.h"
#include "Params.h"
#include "Cell.h"


///
/// Rotation
///

//- rotate coordinate system clockwise
numvector<double, PhysDim> rotate(const numvector<double, PhysDim>& sol, const Point& n);

//- rotate coordinate system counter-clockwise
numvector<double, PhysDim> inverseRotate(const numvector<double, PhysDim>& sol, const Point& n);


///
/// Integration
///

//- 2D Gauss integration of scalar function
double integrate( const Cell& cell, const std::function<double(const Point &)>& f);

//- 2D Gauss integration of vector function
numvector<double, dimPh> integrate(const Cell& cell, const std::function<numvector<double, dimPh>(const Point&)>& f);

///
/// MPI operations
///

int localNumber(std::vector<int> globalNumbers, int curNum);


#endif
