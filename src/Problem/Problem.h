#ifndef PROBLEM_H
#define PROBLEM_H

#include <math.h>
#include <functional>
#include "numvector.h"
#include "Point.h"

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

    //- Coeffs
    std::vector<numvector<double, 5 * nShapes>> alpha;

public:

    //- Default constructor
    Problem();

    //- Destructor
    ~Problem();

    //- Get actual coeffs
    void getAlpha(std::vector<numvector<double, 5 * nShapes> >& a);

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


