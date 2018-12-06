#ifndef PROBLEM_H
#define PROBLEM_H

#include <math.h>
#include <functional>
#include "numvector.h"
#include "Point.h"

#include "Mesh.h"
//#include "Boundary.h"
#include "TimeClass.h"
#include "defs.h"
#include "Params.h"
#include "Physics.h"

enum CaseInit {
    SodX, SodY, SodDiag, SodCircle, Woodward, BlastCircle,
    Noh, OneTwoThree, forwardStep, doubleMach, acousticPulse,
    monopole, dipole,
};

enum CaseBound { ZeroGradient, Inlet, Outlet};

class Patch;

class Problem
{

public:

    //- Function for initial conditions
    std::function<numvector<double, dimPh>(const Point& r)> init;

    //- Parameters on infinity
    numvector<double, dimPh> infty;

    //- Link to coeffs
    //const std::vector<numvector<double, dimPh * nShapes>>& alpha;

    //- Reference to time and space
    const Time& time;
	const Mesh& M;
    Physics& physics;

    //- Constructor
    Problem (CaseInit task, const Time& t, Physics& ph);
    //Problem (CaseInit task);

    //- Destructor
    ~Problem();

    //- Set initial conditions as functions
    void setInitialConditions(CaseInit task);
	//- Set boundary conditions
    void setBoundaryConditions(CaseBound task, const std::vector<Patch> &patches);

};// end Problem

#endif // PROBLEM_H


