#ifndef PROBLEM_H
#define PROBLEM_H

#include <math.h>
#include <functional>
#include "numvector.h"
#include "Point.h"
#include "BoundarySlip.h"
#include "BoundaryOpen.h"
#include "BoundarySine.h"
#include "BoundarySineDir.h"
#include "BoundaryConstant.h"
#include "defs.h"

//- Square
template<class T>
inline T sqr(T x) {return x*x;}

class Patch;


class Problem
{

public:

    //- Heat capacity ratio
    double cpcv;

    //- Function for initial conditions
    std::function<numvector<double, 5>(const Point& r)> init;

    //- Parameters on infinity
    numvector<double, 5> infty;

    //- Coeffs
    std::vector<numvector<double, 5 * nShapes>> alpha;
    //const std::vector<numvector<double, 5 * nShapes>>& alpha;

    //- Reference to time
    const Time& time;

public:

    //- Constructor
    Problem (std::string caseName, const Time& t);

    //- Destructor
    ~Problem();

    //- Set actual coeffs
    void setAlpha(const std::vector<numvector<double, 5 * nShapes> >& a);

    //- Set initial conditions as functions
    void setInitialConditions(std::string caseName);

    void setBoundaryConditions(std::string caseName, const std::vector<Patch> &patches);

    //- Calculate pressure using conservative variables
    double getPressure(const numvector<double,5>& sol) const;

    //- Compute sound speed inside cell
    double c(const numvector<double, 5>& sol) const;

    //- Compute semisum-averaged sound speed on edge
    double c_av(const numvector<double, 5>& solOne, const numvector<double, 5>& solTwo) const;

    //- Eigenvalues for X direction ( Roe )
    numvector<double, 5> lambdaF_Roe(const numvector<double, 5>& solOne, const numvector<double, 5>& solTwo) const;

    //- Eigenvalues for X direction ( Einfeldt )
    numvector<double, 5> lambdaF_Einfeldt(const numvector<double, 5>& solOne, const numvector<double, 5>& solTwo) const;

    //- Eigenvalues for X direction ( Toro, pressure-based )
    numvector<double, 5> lambdaF_Toro(const numvector<double, 5>& solOne, const numvector<double, 5>& solTwo) const;

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

    //- Left  eigenvectors
    std::pair<numvector<numvector<double, 5>, 5>, numvector<numvector<double, 5>, 5>> getL(const numvector<double, 5>& sol) const;
    numvector<numvector<double, 5>, 5> getL(const numvector<double, 5>& sol, const Point& n) const;

    //- Right eigenvectors
    std::pair<numvector<numvector<double, 5>, 5>, numvector<numvector<double, 5>, 5>> getR(const numvector<double, 5>& sol) const;
    numvector<numvector<double, 5>, 5> getR(const numvector<double, 5>& sol, const Point& n) const;

};// end Problem

#endif // PROBLEM_H


