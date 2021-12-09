/// ---------------------
/// SolverHybrid class
/// ---------------------
/// Define main function of solution process
///


#ifndef SOLVERHYBRID_H
#define SOLVERHYBRID_H

#include <fstream>
#include "Point.h"

#include "Basis.h"
#include "Solution.h"
#include "Solver.h"
#include "compService.h"
#include "Problem.h"
#include "Flux.h"

/// Proc rank 
extern int myRank;

/// Size of ...
extern int numProcsTotal;

/// Status
extern MPI_Status status;

/// Debug
extern bool debug;

/// Log file to save data
extern std::ofstream logger;

///
/// Computer of right-hand side of equation system
///

class SolverHybrid : public Solver
{

private: 

    /// Constant reference to Riemann solver
    const Flux& flux1;

    /// Constant reference to Riemann solver
    const Flux& flux2;

    /// Hybrid flux
    void getHybridFlux( 
        const numvector<double, dimPh>& solLeft, 
        const numvector<double, dimPh>& solRight,
        const Point& eNormal,
        numvector<double, dimPh>& flux
    );

public:

    /// Constructor
    SolverHybrid( Basis& B, Mesh& msh, Solution& sln,
            Problem &prb, Physics& phs, Flux& flx1, Flux& flx2, Buffers& buf);

    /// Destructor
    ~SolverHybrid() {}

    //------------ RKDG functions

    /// Assemble right-hand side
    virtual std::vector<numvector<double, dimS>> assembleRHS(const std::vector<numvector<double, dimS>>& SOL);

};

#endif // SolverHybrid_H
