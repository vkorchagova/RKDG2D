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

    //- Coeffs for previous time step
    std::vector<numvector<double, 5 * nShapes>> alphaPrev;

    //- Coeffs for current time step
    std::vector<numvector<double, 5 * nShapes>> alphaNext;

public:

    //- Default constructor
    Solver() {}

    //- Construct using mesh and problem
    Solver(Mesh2D& msh,Problem &prb) { mesh = &msh; problem = &prb; alphaPrev.resize(mesh->nCells); alphaNext.resize(mesh->nCells); }

    // Copy constructor
    Solver(const Solver& rhs) { mesh = rhs.mesh; problem = rhs.problem; alphaPrev = rhs.alphaPrev; alphaNext = rhs.alphaNext; }

    //- Overloaded "=" operator
    Solver& operator=(const Solver& rhs) { mesh = rhs.mesh; problem = rhs.problem; alphaPrev = rhs.alphaPrev; alphaNext = rhs.alphaNext; return *this; }

    //- Destructor
    ~Solver() {}

    //// RKDG methods

    //- Set functions for boundary conditions
    void initBoundaryConditions();

    //- Calculate coeffs with initial conditions
    void setInitialConditions();

    //- Set flux
    void initFluxes(Flux& flux);

    //- Assemble right hand side
    void assembleRHS( std::vector<numvector<double, 5 * nShapes>>& alpha);

    //// Other methods

    //- Output for coeffs
    void write(std::ostream& writer, const numvector<double,5*nShapes>& coeffs);
};

#endif // SOLVER_H
