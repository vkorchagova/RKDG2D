#include <vector>
#include "numvector.h"
#include "defs.h"
#include "Solver.h"
#include "Limiter.h"
#include "TimeClass.h"

///
/// Abstract class for explicit methods of ODE system solution
///

class DDT
{

protected:

    /// Reference to solver
    Solver& solver;

    /// Reference to limiter
    Limiter& limiter;

    /// Reference to time
    Time& time;

    /// Order of accuracy
    int order;

public:

    /// Constructor
    DDT(int o, Solver& s, Limiter& l, Time& t) : order(o), solver(s), limiter(l), time(t) {};

    /// Update time step
    ///
    /// @param  yOld    solution coefficients in all mesh cells from previous time point
    /// @param  tau     time step
    ///
    /// @return solution coefficients in all mesh cells for current time point 
    virtual std::vector<numvector<double, 5*nShapes>> update\
        (std::vector<numvector<double, 5*nShapes>> yOld, double tau) = 0;

};
