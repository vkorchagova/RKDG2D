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

//- Proc rank 
extern int myRank;

//- Size of ...
extern int numProcsTotal;

//- Status
extern MPI_Status status;

//- Request
extern MPI_Request request;


class Solver
{

public:

	/// References
	//- Reference to the basis
	const Basis& B;

	//- Reference to teh mesh
	const Mesh& M;
	
	//- Reference to teh solution
	Solution& sln;

    //- Reference to problem
    const Problem& prb;

    //- Reference to physics
    const Physics& phs;

    //- Reference to flux
    const Flux& flux;

    //- Reference to MPI buffers
    Buffers& buf;

    std::vector<Boundary> bc;

	/// Variables
	//- Max speed buffer for the Courant condition
	double MaxSpeed;

public:

    //- Constructor
    Solver( Basis& B, Mesh& msh, Solution& sln,
			Problem &prb, Physics& phs, Flux& flx, Buffers& buf);

    //- Destructor
    ~Solver() {}

    //// RKDG methods

    //- Calculate coeffs with initial conditions
    void setInitialConditions();

    //- Run case from define set of coefficients
    void setDefinedCoefficients(std::string fileName);

	//- Get coefficients of projection of function foo onto cell basis
	numvector<double, dimS> projection(const std::function<numvector<double, dimPh>(const Point& point)>& init, int iCell) const;

    //- Assemble right-hand side
    std::vector<numvector<double, dimS>> assembleRHS(const std::vector<numvector<double, dimS>>& SOL);

    //- Correct alpha coeffs in case of non-orthogonal basis functions
	numvector<double, dimS> correctNonOrthoCell(const numvector<double, dimS>& rhs, const std::vector<std::vector<double>>& gramian) const;
	std::vector<numvector<double, dimS>> correctNonOrtho(const std::vector<numvector<double, dimS>>& alpha) const;

	//- Reconstruct SLAE RHS after limitation in case of non-orthogonal functions :TODO choose the necessary version
	numvector<double, dimS> correctPrevIterCell(const numvector<double, dimS>& alphaCorr, const std::vector<std::vector<double>>& gramian) const;
	std::vector<numvector<double, dimS>> correctPrevIter(const std::vector<numvector<double, dimS>>& alpha) const;

    /// MPI functions

    //- Exchange data
    void dataExchange();

    //- Collect solution
    void collectSolution();


};

#endif // SOLVER_H
