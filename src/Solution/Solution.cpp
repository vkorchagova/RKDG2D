#include "Solution.h"

#include <iostream>

using namespace std;

///

Solution::Solution(Basis& Bas) : B(Bas)
{

}

///

numvector<double, dimPh> Solution::reconstruct(int iCell, const Point& point ) const
{
    numvector<double, dimPh> sol(0.0);
    vector<double> phip(nShapes);

////#pragma omp simd
    for (int j = 0; j < nShapes; ++j)
        phip[j] = B.phi[j](iCell, point);

    for (int i = 0; i < dimPh; ++i)
    {
        for (int j = 0; j < nShapes; ++j)
            sol[i] += phip[j] * SOL[iCell][i * nShapes + j];
    }

    return sol;

} // end reconstruct
numvector<double, dimPh> Solution::reconstruct(int iCell, const Point& point, const numvector<double, dimS>& SOL) const
{
	numvector<double, dimPh> sol(0.0);
	vector<double> phip(nShapes);

////#pragma omp simd
	for (int j = 0; j < nShapes; ++j)
		phip[j] = B.phi[j](iCell, point);

	for (int i = 0; i < dimPh; ++i)
	{

//#pragma omp simd
		for (int j = 0; j < nShapes; ++j)
			sol[i] += phip[j] * SOL[i * nShapes + j];
	}

	return sol;

} // end reconstruct

double Solution::reconstruct(int iCell, const Point& point, Variables var ) const
{
    double sol(0.0);

////#pragma omp simd    
    for (int j = 0; j < nShapes; ++j)
        sol += B.phi[j](iCell, point) * SOL[iCell][var * nShapes + j];

    return sol;
} // end reconstruct
double Solution::reconstruct(int iCell, const Point& point, const numvector<double, dimS>& SOL, Variables var ) const
{
    double sol(0.0);

////#pragma omp simd
    for (int j = 0; j < nShapes; ++j)
        sol += B.phi[j](iCell, point) * SOL[var * nShapes + j];

    return sol;
} // end reconstruct
