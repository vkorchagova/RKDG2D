#ifndef PROBLEM_H
#define PROBLEM_H

#include <math.h>
#include <functional>
#include "numvector.h"
#include "Point.h"
#include "Mesh.h"
#include "TimeControl.h"
#include "defs.h"
#include "Params.h"

#include "BoundarySlip.h"
#include "BoundaryOpen.h"
#include "BoundaryConstant.h"


class Patch;

class Problem
{

public:

    /// Heat capacity defined by user
    double cpcv;


    /// Function for initial conditions
    std::function<numvector<double, dimPh>(const Point& r)> init;

    /// Vector of boundary conditions
    std::vector<std::shared_ptr<Boundary>> bc;

    /// Parameters on infinity
    numvector<double, dimPh> infty;

    /// Link to coeffs
    //const std::vector<numvector<double, dimPh * nShapes>>& alpha;

    /// Reference to time and space
    const TimeControl& T;
	const Mesh& M;

    /// Constructor
    Problem (CaseInit task, const Mesh& m, const TimeControl& t);

    /// Destructor
    ~Problem();

    /// Set initial conditions as functions
    void setInitialConditions(CaseInit task);

	/// Set boundary conditions
    void setBoundaryConditions(CaseInit task);

};// end Problem

#endif // PROBLEM_H


