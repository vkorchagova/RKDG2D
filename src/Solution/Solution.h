#ifndef SOLUTION_H
#define SOLUTION_H

#include <math.h>
#include <functional>
#include <memory>
#include "numvector.h"

#include "defs.h"
#include "Basis.h"
#include "Physics.h"



class Solution
{

public:   
	
    //- The very Coeffs
    std::vector<numvector<double, dimS>> SOL;
	//- And its auxilliary copy for any technical needs
	std::vector<numvector<double, dimS>> SOL_aux;

public:

    //- Constructor
	Solution();

    //- Destructor
    ~Solution();

	//- Reconstruct solution via inner coeffs
	numvector<double, dimPh> reconstructSolution(const Point& point) const;
	double reconstructSolution(const Point& point, Variables var) const;

    /*numvector<double, dimPh> reconstructSolution(const std::shared_ptr<Point> point) const
	{
		return reconstructSolution(*point);
	}
	double reconstructSolution(const std::shared_ptr<Point> point, Variables var) const
	{
		return reconstructSolution(*point, var);
	}*/

	//- Reconstruct solution using given coeffs
	numvector<double, dimPh> reconstructSolution(const Point& point, const numvector<double, dimS>& SOL) const;
	double reconstructSolution(const Point& point, const numvector<double, dimS>& SOL, Variables var) const;
	/*numvector<double, dimPh> reconstructSolution(const std::shared_ptr<Point> point, const numvector<double, dimS>& SOL) const
	{
		return reconstructSolution(*point, SOL);
	}
	double reconstructSolution(const std::shared_ptr<Point> point, const numvector<double, dimS>& SOL, Variables var) const
	{
		return reconstructSolution(*point, SOL, var);
	}*/

	///
	/// Other methods
	///

	//- Output for coeffs
	void write(std::string fileName, const std::vector<numvector<double, dimS>>& coeffs) const;

	//- Output for VTK solution
	void writeSolutionVTK(std::string fileName) const;

};// end Solution

#endif // SOLUTION_H


