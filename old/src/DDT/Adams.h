/// 
/// Explicit Adams methods
/// 

#include "RungeKutta.h"
#include <deque>

class Adams : public RungeKutta
{

    //- history of rhs
    std::deque<std::vector<numvector<double, 5 * nShapes>>> rhsHistory;

    //- Adams coeffs
    std::vector<double> b;
    
    //- Set coeffs for defined order
    void setAdamsCoeffs();
    
public:

    //- Constructor
    Adams(int o, Solver& s, Limiter& l, Time& t);
    
    //- Destructor
    ~Adams() {};

    //- Compute time step
    virtual std::vector<numvector<double, 5*nShapes>> update\
         (std::vector<numvector<double, 5*nShapes>> yOld, double tau);

};