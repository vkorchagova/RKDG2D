#ifndef COMPSERVICE_H_
#define COMPSERVICE_H_

#include "numvector.h"
#include "Params.h"
#include "Point.h"
#include "Cell.h"
#include "Patch.h"

///
/// Special functions for computations
///

///
/// Rotation
///

/// rotate coordinate system clockwise
numvector<double, dimPh> rotate(const numvector<double, dimPh>& sol, const Point& n);
Point rotate(const Point& v, const Point& n);

/// rotate coordinate system counter-clockwise
numvector<double, dimPh> inverseRotate(const numvector<double, dimPh>& sol, const Point& n);
Point inverseRotate(const Point& v, const Point& n);


///
/// Integration
///

/// 1D Gauss integration of scalar function
double integrate( const Edge& edge, const std::function<double(const Point &)>& f);

/// 1D Gauss integration of scalar function by defined values in Gauss points
double integrate( const Edge& edge, const std::vector<double>& f);

/// 1D Gauss integration of vector function
numvector<double, dimPh> integrate(const Edge& edge, const std::function<numvector<double, dimPh>(const Point&)>& f);

/// 1D Gauss integration of vector function by defined values in Gauss points
numvector<double, dimPh> integrate(const Edge& edge, const std::vector<numvector<double, dimPh>>& f);


/// 2D Gauss integration of scalar function
double integrate( const Cell& cell, const std::function<double(const Point &)>& f);

/// 2D Gauss integration of vector function
numvector<double, dimPh> integrate(const Cell& cell, const std::function<numvector<double, dimPh>(const Point&)>& f);


///
/// MPI operations
///

/// Get local cell number on proc according to global cell number
int localNumber(std::vector<int>& globalNumbers, int curNum);

int getPatchByName(std::vector<Patch>& patches, std::string& pName);
int getPatchByName(std::vector<ProcPatch>& patches, std::string& pName);


///
/// Viscosity
///

numvector<double, dimGrad> outerProductArtificial(const numvector<double, dimPh>& sol, const Point& n);


#endif
