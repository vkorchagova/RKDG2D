/// ---------------------
/// Solver class
/// ---------------------
/// Define main function of solution process
///


#ifndef SOLVER_H
#define SOLVER_H

#include <fstream>
#include "Point.h"

#include "Basis.h"
#include "Solution.h"
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

class Solver
{
    /// RHS
    std::vector<numvector<double, dimS>> rhs; // the same length as SOL
    
    /// List of numerical fluxes in Gauss points
    std::vector<std::vector<numvector<double, dimPh>>> numFluxes;

public:

    /// Constant reference to the basis
    const Basis& B;

    /// Constant reference to teh mesh
    const Mesh& M;
    
    /// Reference to teh solution
    Solution& sln;

    /// Constant reference to problem
    const Problem& prb;

    /// Constant reference to physics
    const Physics& phs;

    /// Constant reference to Riemann solver
    const Flux& flux;

    /// Reference to MPI buffers
    Buffers& buf;

    /// Vector of boundary conditions
    //std::vector<Boundary> bc;

    /// Max speed buffer for the Courant condition (???)
    double MaxSpeed;

public:

    /// Constructor
    Solver( Basis& B, Mesh& msh, Solution& sln,
            Problem &prb, Physics& phs, Flux& flx, Buffers& buf);

    /// Destructor
    ~Solver() {}

    //------------ RKDG functions

    /// Calculate coeffs with initial conditions
    void setInitialConditions();

    /// Run case from define set of coefficients
    void restart(std::string fileName);

    /// Assemble right-hand side
    std::vector<numvector<double, dimS>> assembleRHS(const std::vector<numvector<double, dimS>>& SOL);

    /// Correct alpha coeffs in case of non-orthogonal basis functions
    numvector<double, dimS> correctNonOrthoCell(const numvector<double, dimS>& rhs, const std::vector<std::vector<double>>& gramian) const;
    std::vector<numvector<double, dimS>> correctNonOrtho(const std::vector<numvector<double, dimS>>& alpha) const;

    /// Reconstruct SLAE RHS after limitation in case of non-orthogonal functions :TODO choose the necessary version
    numvector<double, dimS> correctPrevIterCell(const numvector<double, dimS>& alphaCorr, const std::vector<std::vector<double>>& gramian) const;
    std::vector<numvector<double, dimS>> correctPrevIter(const std::vector<numvector<double, dimS>>& alpha) const;

    //------------ MPI functions

    /// Exchange data
    void dataExchange();

    /// Collect solution
    void collectSolution();
    void collectSolutionForExport();


};

#endif // SOLVER_H
