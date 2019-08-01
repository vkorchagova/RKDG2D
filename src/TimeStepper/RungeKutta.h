#include "TimeStepper.h"

/// 
/// Explicit Runge --- Kutta numerical scheme
/// 

class RungeKutta : public TimeStepper
{
    /// Alpha Butcher coeffs for RK studies
    std::vector<double> alpha;

    /// Beta Butcher coeffs for RK studies
    std::vector<std::vector<double>> beta;

    /// set Butcher coeffs for given order
    void setButcherTable();

public:

    /// Constructor
    RungeKutta(int o,  Basis& b, Solver& s, Solution& ss, Limiter& l, TimeControl& t);
    
    /// Destructor
    ~RungeKutta() {};

    /// Compute time step
    virtual void Tstep();

};