// Time derivative implementation

#include <vector>
#include <deque>
#include "Solver.h"
#include "Boundary.h"
#include "Limiter.h"
#include "TimeControl.h"

/// Debug
extern bool debug;

/// Log file to save data
extern std::ofstream logger;

///
/// Abstract class for time step numerical scheme
///

class TimeStepper
{

protected:

    /// Reference to solver
    Solver& slv;

    /// Reference to teh solution
    Solution& sln;

    /// Reference to basis
    Basis& basis;

    /// Reference to the Full Pack of Boundary Conditions
    std::vector<std::shared_ptr<Boundary>>& bc;

    /// Reference to limiter
    Limiter& lmt;

    /// Reference to time
    TimeControl& T;

    /// Order of accuracy
    int order;

    /// Number of RK/Adams steps --- TODO:  is once needed to separate from order!
    int nStages;

public:

    /// Constructor
    TimeStepper(int o, Basis& b, Solver& s, Solution& ss, std::vector<std::shared_ptr<Boundary>>& bond, Limiter& l, TimeControl& t) : order(o), basis(b), slv(s), sln(ss), bc(bond), lmt(l), T(t) {};

    /// update time step
    virtual void Tstep() = 0;

};
