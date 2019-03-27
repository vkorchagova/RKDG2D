/// 
/// Explicit Runge --- Kutta methods
/// 

#include "TimeStepper.h"

class RungeKutta : public TimeStepper
{
    

    /// Butcher coeffs for RK studies
    std::vector<double> alpha;
    std::vector<std::vector<double>> beta;

    /// set Butcher coeffs for given order
    void setButcherTable();

public:

    /// Constructor
    RungeKutta(int o,  Basis& b, Solver& s, Solution& ss, std::vector<std::shared_ptr<Boundary>>& bond, Limiter& l, TimeControl& t);
    
    /// Destructor
    ~RungeKutta() {};

    /// Compute time step
    virtual void Tstep();

};