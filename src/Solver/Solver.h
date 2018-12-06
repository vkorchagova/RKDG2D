/// ---------------------
/// Solver class
/// ---------------------
/// Define main function of solution process
///


#ifndef SOLVER_H
#define SOLVER_H

#include <fstream>

#include "Basis.h"
#include "compService.h"
#include "Problem.h"
#include "Flux.h"
#include "Limiter.h"

class Solver
{

public:

	//- Reference to the basis
	const Basis& B;

	//- Reference to teh mesh
	const Mesh& M;
	
	//- Reference to teh solution
	Solution& sln;

    //- Reference to Gas problem
    const Problem& prb;

    //- Reference to flux
    const Flux& flux;

public:

    //- Constructor
    Solver( Basis& B, Mesh& msh, Solution& sln,
			Problem &prb, Flux& flx);

    //- Destructor
    ~Solver() {}

    //// RKDG methods

    //- Calculate coeffs with initial conditions
    void setInitialConditions();

    //- Run case from define set of coefficients
    void setDefinedCoefficients(std::string fileName);

	//- Get coefficients of projection of function foo onto cell basis
	numvector<double, dimS> projection(const std::function<numvector<double, dimPh>(const Point& point)>& init, const Cell& cell) const;

    //- Assemble right-hand side
    std::vector<numvector<double, dimS>> assembleRHS(const std::vector<numvector<double, dimS>>& SOL);

    //- Correct alpha coeffs in case of non-orthogonal basis functions
	numvector<double, dimS> correctNonOrtho(const numvector<double, dimS>& alpha, const Cell& cell) const;

	//- Reconstruct SLAE RHS after limitation in case of non-orthogonal functions
	numvector<double, dimS> correctPrevIter(const numvector<double, dimS>& rhs, const Cell& cell) const;

};

#endif // SOLVER_H
