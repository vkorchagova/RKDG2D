#ifndef COMPSERVICE_H_
#define COMPSERVICE_H_

#include "numvector.h"
#include "Params.h"
#include "Point.h"
#include "Cell.h"
#include "Patch.h"


///
/// Rotation
///

//- rotate coordinate system clockwise
numvector<double, PhysDim> rotate(const numvector<double, PhysDim>& sol, const Point& n);
Point rotate(const Point& v, const Point& n);

//- rotate coordinate system counter-clockwise
numvector<double, PhysDim> inverseRotate(const numvector<double, PhysDim>& sol, const Point& n);
Point inverseRotate(const Point& v, const Point& n);


///
/// Integration
///

//- 1D Gauss integration of scalar function
double integrate( const Edge& edge, const std::function<double(const Point &)>& f);

//- 1D Gauss integration of scalar function by defined values in Gauss points
double integrate( const Edge& edge, const std::vector<double>& f);

//- 1D Gauss integration of vector function
numvector<double, PhysDim> integrate(const Edge& edge, const std::function<numvector<double, PhysDim>(const Point&)>& f);

//- 1D Gauss integration of vector function by defined values in Gauss points
numvector<double, PhysDim> integrate(const Edge& edge, const std::vector<numvector<double, PhysDim>>& f);


//- 2D Gauss integration of scalar function
double integrate( const Cell& cell, const std::function<double(const Point &)>& f);

//- 2D Gauss integration of vector function
numvector<double, PhysDim> integrate(const Cell& cell, const std::function<numvector<double, PhysDim>(const Point&)>& f);


///
/// MPI operations
///

int localNumber(std::vector<int>& globalNumbers, int curNum);

int getPatchByName(std::vector<Patch>& patches, std::string& pName);
int getPatchByName(std::vector<ProcPatch>& patches, std::string& pName);


#endif
