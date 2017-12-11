/// ---------------------
/// Solver class
/// ---------------------
/// Define main function of solution process
///


#ifndef SOLVER_H
#define SOLVER_H

#include <fstream>

#include "Mesh2D.h"
#include "Problem.h"

class Solver
{

public:

    //- Pointer to 2D mesh
    Mesh2D* mesh;

    //- Pointer to Gas2D problem
    Problem* problem;

public:

    //- Default constructor
    Solver() {}

    //- Construct using mesh and problem
    Solver(Mesh2D& msh,Problem &prb) { mesh = &msh; problem = &prb; }

    // Copy constructor
    Solver(const Solver& rhs) { mesh = rhs.mesh; problem = rhs.problem; }

    //- Overloaded "=" operator
    Solver& operator=(const Solver& rhs) { mesh = rhs.mesh; problem = rhs.problem; return *this; }

    //- Destructor
    ~Solver() {}

    //// RKDG methods

    //- Set initial conditions
    void setInitialConditions();

    //- Assemble right hand side
    void assembleRHS();

    //// Other methods

    //- Output for coeffs
    void write(std::ostream& writer, const numvector<double,5*nShapes>& coeffs);
};

#endif // SOLVER_H
