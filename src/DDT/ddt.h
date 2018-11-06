// Time derivative implementation

#include <vector>
#include "numvector.h"
#include "defs.h"
#include "Solver.h"
#include "Limiter.h"
#include "TimeClass.h"

class DDT
{

protected:

    //- Reference to solver
    Solver& solver;

    //- Reference to limiter
    Limiter& limiter;

    //- Reference to time
    Time& time;

    //- Order of accuracy
    int order;

    //- Number of RK steps
    int nSteps;

public:

    //- Constructor
    DDT(int o, Solver& s, Limiter& l, Time& t) : order(o), solver(s), limiter(l), time(t) {};

    //- update time step
    virtual std::vector<numvector<double, 5*nShapes>> update\
        (std::vector<numvector<double, 5*nShapes>> yOld, double tau) = 0;

};
