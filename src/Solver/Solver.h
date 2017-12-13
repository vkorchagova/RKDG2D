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
    Solver( Mesh2D& msh, Problem &prb) { mesh = &msh; problem = &prb; alphaPrev.resize(mesh->nCells); alphaNext.resize(mesh->nCells); }

    //- Destructor
    ~Solver() {}

    //// RKDG methods

    //- Set functions for boundary conditions
    void initBoundaryConditions() const;

    //- Calculate coeffs with initial conditions
    void setInitialConditions() const;

    //- Set flux
    void initFluxes(const Flux& flux) const;

    //- Assemble right hand side
    void assembleRHS( const std::vector<numvector<double, 5 * nShapes>>& alpha) const;

    //// Other methods

    //- Output for coeffs
    void write(std::ostream& writer, const numvector<double,5*nShapes>& coeffs) const;
};

#endif // SOLVER_H
