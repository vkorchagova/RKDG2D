#ifndef SOLUTION_H
#define SOLUTION_H

#include <math.h>
#include <functional>
#include <memory>
#include "numvector.h"

#include "defs.h"
#include "Basis.h"
#include "Params.h"


class Solution
{

public:   
	
    //- The very Coeffs
    std::vector<numvector<double, dimS>> SOL;

	//- And its auxilliary copy for any technical needs
	std::vector<numvector<double, dimS>> SOL_aux;

	//- Reference to the basis
	const Basis& B;

public:

    //- Constructor
	Solution(Basis& bas);

    //- Destructor
	~Solution() {}

	//- Reconstruct solution at the point
	numvector<double, dimPh> reconstruct(int iCell, const Point& point) const;
	double reconstruct(int iCell, const Point& point, Variables var) const;
	/*numvector<double, dimPh> reconstruct(const std::shared_ptr<Point> point) const
	{
		return reconstruct(*point);
	}
	double reconstruct(const std::shared_ptr<Point> point, Variables var) const
	{
		return reconstruct(*point, var);
	}*/

	//- Reconstruct solution using given coeffs
	numvector<double, dimPh> reconstruct(int iCell, const Point& point, const numvector<double, dimS>& SOL) const;
	double reconstruct(int iCell, const Point& point, const numvector<double, dimS>& SOL, Variables var) const;
	/*numvector<double, dimPh> reconstruct(const std::shared_ptr<Point> point, const numvector<double, dimS>& SOL) const
	{
		return reconstruct(*point, SOL);
	}
	double reconstruct(const std::shared_ptr<Point> point, const numvector<double, dimS>& SOL, Variables var) const
	{
		return reconstruct(*point, SOL, var);
	}*/

	///
	/// Other methods
	///

	//- Output for coeffs
	//void write(std::string fileName, const std::vector<numvector<double, dimS>>& coeffs) const;

	//- Output for VTK solution
	//void writeSolutionVTK(std::string fileName) const;

};// end Solution

#endif // SOLUTION_H


