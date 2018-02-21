#ifndef PROBLEM_H
#define PROBLEM_H

#include <math.h>
#include <functional>
#include "numvector.h"
#include "Point.h"
#include "BoundarySlip.h"
#include "defs.h"

//- Square
template<class T>
inline T sqr(T x) {return x*x;}

class Patch;


class Problem
{

public:

    //- Heat capacity ratio
    const double cpcv = 1.4;

    //- Function for initial conditions
    std::function<numvector<double, 5>(const Point& r)> init;

    //- Parameters on infinity
    numvector<double, 5> infty;

    //- Coeffs
    std::vector<numvector<double, 5 * nShapes>> alpha;
    //const std::vector<numvector<double, 5 * nShapes>>& alpha;

public:

    //- Default constructor
    Problem ();

    //- Destructor
    ~Problem();

    //- Set actual coeffs
    void setAlpha(const std::vector<numvector<double, 5 * nShapes> >& a);

    //- Set initial conditions as functions
    void setInitialConditions();

    void setBoundaryConditions(const std::vector<Patch> &patches);

    //- Calculate pressure using conservative variables
    double getPressure(const numvector<double,5>& sol) const;

    //- Compute sound speed inside cell
    double c(const numvector<double, 5>& sol) const;

    //- Compute semisum-averaged sound speed on edge
    double c_av(const numvector<double, 5>& solOne, const numvector<double, 5>& solTwo) const;

    //- Eigenvalues for X direction ( Roe )
    numvector<double, 5> lambdaF_Roe(const numvector<double, 5>& solOne, const numvector<double, 5>& solTwo) const;

    //- Eigenvalues for X direction ( semisum )
    numvector<double, 5> lambdaF_semisum(const numvector<double, 5>& solOne, const numvector<double, 5>& solTwo) const;

    //- Eigenvalues for X direction
    numvector<double, 5> lambdaF(const numvector<double, 5>& solOne, const numvector<double, 5>& solTwo) const;

    //- Eigenvalues for Y direction ( semisum )
    numvector<double, 5> lambdaG(const numvector<double, 5>& solOne, const numvector<double, 5>& solTwo) const;

    //- Calculate fluxes in x direction
    numvector<double, 5> fluxF(const numvector<double, 5>& sol) const;

    //- Calculate fluxes in y direction
    numvector<double, 5> fluxG(const numvector<double, 5>& sol) const;

};// end Problem

#endif // PROBLEM_H


