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
#include "Physics.h"

#include "BoundarySlip.h"
#include "BoundaryOpen.h"
#include "BoundaryConstant.h"
#include "BoundaryNonSlip.h"
#include "BoundaryInlet.h"

class Patch;

///
/// Store of mesh, time, initial and boundary conditions ("world") 
///

class Problem
{

private:

    /// Constant reference to time 
    const TimeControl& T;

    /// Constant reference to mesh
    const Mesh& M;

    /// Reference to physics
    Physics& phs;

    /// Set initial conditions as functions
    void setInitialConditions(CaseInit task);

    /// Set boundary conditions
    void setBoundaryConditions(CaseInit task);


public:

    /// Function for initial conditions
    std::function<numvector<double, dimPh>(const Point& r)> init;

    /// Function for explicit source
    std::function<numvector<double, dimPh>(const numvector<double, dimPh> sol, const Point& r)> source;

    /// Vector of boundary conditions
    std::vector<std::shared_ptr<Boundary>> bc;

    /// Constructor
    Problem (CaseInit task, const Mesh& m, const TimeControl& t, Physics& phs);

    /// Destructor
    ~Problem();

};// end Problem

#endif // PROBLEM_H


