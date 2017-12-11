#ifndef PROBLEM_H
#define PROBLEM_H

#include "Mesh2D.h"
#include "gaussintegrator.h"
//#include "fluxllf.h"

#include <math.h>
#include <functional>

class Mesh2D;

//- Number of basis functions
const int nShapes = 3;

class Problem
{

public:


    //- Heat capacity ratio
    const double cpcv = 1.4;

    //- Function for initial conditions
    std::function<numvector<double, 5>(const Point& r)> init;

    //- Parameters on infinity
    numvector<double,5> infty;

public:

    //- Default constructor
    Problem();

    //- Destructor
    ~Problem();

    //- Calculate pressure using conservative variables
    double getPressure(numvector<double,5> sol);

    //- Compute sound speed inside cell
    double c(numvector<double, 5> sol);

    //- Compute averaged sound speed on edge
    double c_av(numvector<double, 5> solOne, numvector<double, 5> solTwo);

    //- Eigenvalues for X direction ( estimation!!! )
    numvector<double, 5> lambdaF(numvector<double, 5> solOne, numvector<double, 5> solTwo);

    //- Eigenvalues for Y direction ( estimation!!! )
    numvector<double, 5> lambdaG(numvector<double, 5> solOne, numvector<double, 5> solTwo);
	
//    //- Reconstruct solution by coeffs and basis functions
//    numvector<double, 5> reconstructSolution(const numvector<double, \
//                                            nShapes * 5>& alpha, \
//											numvector<double, 2>& coord, \
//											int iCell);


    //- Calculate fluxes in x direction
    numvector<double, 5> fluxF(numvector<double, 5> sol);

    //- Calculate fluxes in y direction
    numvector<double, 5> fluxG(numvector<double, 5> sol);

};// end Problem

#endif // PROBLEM_H


