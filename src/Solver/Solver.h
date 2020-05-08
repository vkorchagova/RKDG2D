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

private: 
    
    /// RHS
    std::vector<numvector<double, dimS>> rhs; // the same length as SOL
    
    /// GradU
    std::vector<numvector<double, dimGradCoeff>> gradU;
    
    /// List of numerical fluxes in Gauss points
    std::vector<std::vector<numvector<double, dimPh>>> numFluxes;

    /// List of numerical fluxes in Gauss points for Grad U
    std::vector<std::vector<numvector<double, dimGrad>>> HnumFluxes;

    /// List of viscous fluxes in Gauss points
    std::vector<std::vector<numvector<double, dimPh>>> VnumFluxes;

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

    /// Constant reference to Riemann solver
    const Flux& vflux;

    /// Reference to MPI buffers
    Buffers& buf;

public:

    /// Constructor
    Solver( Basis& B, Mesh& msh, Solution& sln,
            Problem &prb, Physics& phs, Flux& flx, Flux& vflx, Buffers& buf);

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
    std::vector<numvector<double, dimS>> correctNonOrtho(const std::vector<numvector<double, dimS>>& alpha) const;

    /// Reconstruct SLAE RHS after limitation in case of non-orthogonal functions :TODO choose the necessary version
    std::vector<numvector<double, dimS>> correctPrevIter(const std::vector<numvector<double, dimS>>& alpha) const;

    /// Function to compute grad of conservative variables
    void computeGradU(const std::vector<numvector<double, dimS>>& SOL, std::vector<numvector<double, dimGradCoeff>>& S);

    /// Functon to compute viscous fluxes in Gauss points for RHS
    void computeVnumFluxes();
    //------------ MPI functions

    /// Exchange data
    void dataExchange();

    /// Collect solution
    void collectSolution();
    void collectSolutionForExport();


};

#endif // SOLVER_H
