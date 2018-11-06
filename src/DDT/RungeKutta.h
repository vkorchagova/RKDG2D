/// 
/// Explicit Runge --- Kutta methods
/// 

#include "ddt.h"

class RungeKutta : public DDT
{
    //- rhs for RK studies
    std::vector<std::vector<numvector<double, 5*nShapes>>> k;

    //- Butcher coeffs for RK studies
    std::vector<double> alpha;
    std::vector<std::vector<double>> beta;

    //- set Butcher coeffs for given order
    void setButcherTable();

public:

    //- Constructor
    RungeKutta(int o, Solver& s, Limiter& l, Time& t);
    
    //- Destructor
    ~RungeKutta() {};

    //- Compute time step
    virtual std::vector<numvector<double, 5*nShapes>> update\
         (std::vector<numvector<double, 5*nShapes>> yOld, double tau);

};