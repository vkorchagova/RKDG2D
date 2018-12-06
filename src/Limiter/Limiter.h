#ifndef LIMITER_H
#define LIMITER_H

#include "Indicator.h"
#include "Physics.h"
#include "Solution.h"
#include "Mesh.h"
#include "defs.h"
#include <vector>

class Limiter
{

protected:

    //- Discontinuities checker
    const Indicator& indicator;
    //- Physics
	const Physics& phs;
	//- Solution
	const Solution& sln;
	//- Mesh
	const Mesh& M;

    //- Number of limitation steps
    static const int nIter = 2;

public:

    //- Construct by indicator and problem
    Limiter(const Indicator& ind, Physics& phys, Solution& sol, Mesh& mesh) : indicator(ind), phs(phys), sln(sol), M(mesh) {}

    //- Limit solution gradients
    virtual void limit(std::vector<numvector<double, dimS>>& SOL) = 0;

};

#endif // LIMITER_H
