// Time derivative implementation

#include <vector>
#include <deque>
#include "Solver.h"
#include "Boundary.h"
#include "Limiter.h"
#include "TimeControl.h"

class TimeStepper
{

protected:

	//- Reference to solver
	Solver& slv;

	//- Reference to teh solution
	Solution& sln;

    Basis& basis;

	//- Reference to the Full Pack of Boundary Conditions
	std::vector<std::shared_ptr<Boundary>>& bc;

	//- Reference to limiter
	Limiter& lmt;

    //- Reference to time
    TimeControl& T;

    //- Order of accuracy
    int order;

    //- Number of RK/Adams steps --- TODO:  is once needed to separate from order!
    int nStages;

	//- An array for the accumulative stuff: inner stages in RK and previous RHSs in Adams
	//std::deque<std::vector<numvector<double, dimS>>> Arr;

public:

    //- Constructor
	TimeStepper(int o, Basis& b, Solver& s, Solution& ss, std::vector<std::shared_ptr<Boundary>>& bond, Limiter& l, TimeControl& t) : order(o), basis(b), slv(s), sln(ss), bc(bond), lmt(l), T(t) {};

    //- update time step
    virtual void Tstep() = 0;

};
