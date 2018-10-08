#include "ddt.h"

/// 
/// Explicit Runge --- Kutta methods for solution of ODE systems
///

class RungeKutta : public DDT
{
    /// RHS for all RK studies
    std::vector<std::vector<numvector<double, 5*nShapes>>> k;

    /// Butcher alpha coeffs for RK studies
    std::vector<double> alpha;

    /// Butcher alpha coeffs for RK studies
    std::vector<std::vector<double>> beta;

    /// set Butcher coeffs for given order
    void setButcherTable();

    /// Number of RK steps
    int nSteps;

public:

    /// Constructor
    RungeKutta(int o, Solver& s, Limiter& l, Time& t);
    
    /// Destructor
    ~RungeKutta() {};

    /// Update time step
    ///
    /// @param  yOld    solution coefficients in all mesh cells from previous time point
    /// @param  tau     time step
    ///
    /// @return solution coefficients in all mesh cells for current time point 
    virtual std::vector<numvector<double, 5*nShapes>> update\
         (std::vector<numvector<double, 5*nShapes>> yOld, double tau);

};