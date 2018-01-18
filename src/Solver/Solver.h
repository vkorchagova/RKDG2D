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
#include "Flux.h"

class Solver
{

public:

    //- Reference to 2D mesh
    Mesh2D& mesh;

    //- Reference to Gas2D problem
    Problem& problem;

    //- Reference to flux
    Flux& flux;

    //- Coeffs for previous time step
    std::vector<numvector<double, 5 * nShapes>> alphaPrev;

    //- Coeffs for current time step
    std::vector<numvector<double, 5 * nShapes>> alphaNext;

public:

    //- Construct using mesh and problem
    Solver( Mesh2D& msh, Problem &prb, Flux& flx);

    //- Destructor
    ~Solver() {}

    //// RKDG methods

    //- Set functions for boundary conditions
    void setBoundaryConditions() const;

    //- Calculate coeffs with initial conditions
    void setInitialConditions();

    //- Set mesh pointer in case of DiagProject BC
    void setMeshPointerForDiagBC();

    //- Set flux
    //void initFluxes(const Flux& flux) const;

    //- Assemble right-hand side
    std::vector<numvector<double, 5 * nShapes> > assembleRHS( const std::vector<numvector<double, 5 * nShapes>>& alpha);

    //// Other methods

    //- Output for coeffs
    void write(std::ostream& writer, const std::vector<numvector<double,5*nShapes>>& coeffs) const;
};

#endif // SOLVER_H
