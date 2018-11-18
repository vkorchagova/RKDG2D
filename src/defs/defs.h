#ifndef DEFS_H_
#define DEFS_H_

#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cmath>
#include <functional>

#include "numvector.h"
#include "Params.h"
#include "Point.h"

///
/// Redefinition of the dimensions
///

// Coeffs vector size
const int dimS = PhysDim * nShapes;

// Solution vector size
const int dimPh = PhysDim;


//- rotate coordinate system clockwise
numvector<double, dimPh> rotate(const numvector<double, dimPh>& sol, const Point& n);

//- rotate coordinate system counter-clockwise
numvector<double, dimPh> inverseRotate(const numvector<double, dimPh>& sol, const Point& n);


///
/// Overloading of vector<numvector<double,dimPhys * dimShapes>> 
///

//- Multiplying of the solution onto the number
std::vector<numvector<double, dimS>> operator * (const std::vector<numvector<double, dimS>>& a, const double b);

//- Sum of two vector<numvector>'s
std::vector<numvector<double, dimS>> operator + (const std::vector<numvector<double, dimS>>& b, const std::vector<numvector<double, dimS>>& a);


///
/// Overloading for vector<vector<double>> 
///

// Overload of vvd += vvd
std::vector<std::vector<double>>& operator += (std::vector<std::vector<double>>& a, const std::vector<std::vector<double>>& b);

// Overload of vvd -= vvd
std::vector<std::vector<double>>& operator -= (std::vector<std::vector<double>>& a, const std::vector<std::vector<double>>& b);

// Overload of vvd * double
std::vector<std::vector<double>>& operator *= (std::vector<std::vector<double>>& a, const double b);

// Overload of double * vvd
std::vector<std::vector<double>> operator * (const std::vector<std::vector<double>>& a, const double b);


///
/// Overloading for vector<double> 
///

// Overload of vd *= double
std::vector<double>& operator *= (std::vector<double>& a, const double b);

std::vector<double> operator * (const std::vector<double>& a, const double b);

// Overload of double * vd
std::vector<double> operator * (const double b, const std::vector<double>& a);

// Overload of vd / double
std::vector<double> operator / (const std::vector<double>& a, const double b);

// Overload of vd += vd
std::vector<double>& operator += (std::vector<double>& a, const std::vector<double>& b);

// Overload of vd + vd
std::vector<double> operator + (const std::vector<double>& a, const std::vector<double>& b);

// Overload of vd - vd
std::vector<double> operator - (const std::vector<double>& a, const std::vector<double>& b);

// Overload of vd -= vd
std::vector<double>& operator -= (std::vector<double>& a, const std::vector<double>& b);


#endif
