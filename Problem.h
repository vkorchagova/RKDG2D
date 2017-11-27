#pragma once

#include "Mesh2D.h"
#include "gaussintegrator.h"

#include <math.h>
#include <functional>
#include <fstream>

namespace std
{

class Problem
{

private:


    //- Heat capacity ratio
    const double cpcv = 1.4;

    //- Mesh
    Mesh2D mesh;

    //- Number of basis functions
    static const int nShapes = 3;

    //- File for ofstream
    ofstream writer;

private:

    //- Output for coeffs
    void write(ostream& writer, const numvector<double,5*nShapes>& coeffs);

    //- Calculate pressure using conservative variables
    double getPressure(numvector<double,5> sol);
	
    //- Reconstruct solution by coeffs and basis functions
    numvector<double, 5> reconstructSolution(numvector<double, \
											nShapes + 5>& alpha, \
											numvector<double, 2>& coord, \
											int iCell);


    //- Calculate fluxes in x direction
    numvector<double, 5> fluxF(numvector<double, 5> sol);

    //- Calculate fluxes in y direction
	numvector<double, 5> fluxG(numvector<double, 5> sol);

    

public:

    //- List of basis functions
    function <double(const numvector<double, 2>&, int)> phi[nShapes];

    //- Gradient of basis functions
    function <numvector<double, 2>(const numvector<double, 2>&, int)> gradPhi[3];

	//- Solution coeffs on cells on previous time step
	//  2D array: ncells x (nshapes*5)
	vector<numvector<double, 5 * nShapes>> alphaPrev;

	//- Solution coeffs on cells on next time step
	//  2D array: ncells x (nshapes*5)
    vector<numvector<double, 5 * nShapes>> alphaNext;

	//- Set initial conditions
	void setInitialConditions();


public:
    Problem() {};

    Problem(const Mesh2D& mesh);

    ~Problem();





};// end Problem

} // end namespace std;

