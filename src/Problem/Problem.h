#pragma once

#include "Mesh2D.h"
#include "gaussintegrator.h"
//#include "fluxllf.h"

#include <math.h>
#include <functional>
#include <fstream>



//- Number of basis functions
const int nShapes = 3;

class Mesh2D;

class Problem
{

public:


    //- Heat capacity ratio
    const double cpcv = 1.4;



    //- Mesh
    Mesh2D *mesh;

    //- File for ofstream
    std::ofstream writer;

public:

    //- Default constructor
    Problem() {}

    //- Construct with mesh
    Problem(Mesh2D &mesh);

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
	
    //- Reconstruct solution by coeffs and basis functions
    numvector<double, 5> reconstructSolution(const numvector<double, \
                                            nShapes * 5>& alpha, \
											numvector<double, 2>& coord, \
											int iCell);

    //- Set initial conditions
    void setInitialConditions();


    //- Calculate fluxes in x direction
    numvector<double, 5> fluxF(numvector<double, 5> sol);

    //- Calculate fluxes in y direction
    numvector<double, 5> fluxG(numvector<double, 5> sol);

    //- Output for coeffs
    void write(std::ostream& writer, const numvector<double,5*nShapes>& coeffs);

};// end Problem


