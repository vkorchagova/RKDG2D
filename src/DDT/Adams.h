#include "RungeKutta.h"
#include <deque>

/// 
/// Explicit Adams methods for solution of ODE systems
///

class Adams : public RungeKutta
{

    /// History of rhs
    std::deque<std::vector<numvector<double, 5 * nShapes>>> rhsHistory;

    /// Adams coeffs
    std::vector<double> b;

    /// Number of studies
    int nSteps;
    
    /// Set coeffs for defined order
    void setAdamsCoeffs();
    
public:

    /// Constructor
    Adams(int o, Solver& s, Limiter& l, Time& t);
    
    /// Destructor
    ~Adams() {};

    /// Update time step
    ///
    /// @param  yOld    solution coefficients in all mesh cells from previous time point
    /// @param  tau     time step
    ///
    /// @return solution coefficients in all mesh cells for current time point 
    virtual std::vector<numvector<double, 5*nShapes>> update\
         (std::vector<numvector<double, 5*nShapes>> yOld, double tau);

};