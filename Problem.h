#pragma once

#include "Mesh2D.h"

#include <math.h>
#include<functional>

const int nShapes = 3;


const double gamma = 1.4;

class Problem
{

private:

	double phi0(numvector<double,2> r) { return 0.25; }
	double phi1(numvector<double, 2> r) { return 0.5*sqrt(3)*r[0]; }
	double phi2(numvector<double, 2> r) { return 0.5*sqrt(3)*r[1]; }

	Mesh2D mesh;

private:

	double getPressure(numvector<double,5> sol);
	
	numvector<double, 5> reconstructSolution(numvector<double, \
											nShapes + 5>& alpha, \
											numvector<double, 2>& coord);


	numvector<double, 5> fluxF(numvector<double, 5> sol);
	numvector<double, 5> fluxG(numvector<double, 5> sol);

public:

	//array of basis functions - ?

	//- Solution coeffs on cells on previous time step
	//  2D array: ncells x (nshapes*5)
	vector<numvector<double, 5*nShapes>> alphaPrev;

	//- Solution coeffs on cells on next time step
	//  2D array: ncells x (nshapes*5)
	vector<numvector<double, 5 * nShapes>> alphaNext;

	//- Array of basis functions
    vector<function <double(*)(numvector<double, 2>)>> phi;


public:
	Problem();
	Problem(Mesh2D& mesh);


	~Problem();




};

