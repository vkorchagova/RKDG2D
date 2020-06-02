#ifndef SOLUTION_H
#define SOLUTION_H

#include <math.h>
#include <functional>
#include <memory>

#include "numvector.h"
#include "defs.h"
#include "Basis.h"
#include "Params.h"

///
/// Store of solution as coeffs near basis functions
///

class Solution
{

public:   

    /// Reference to the basis PUBLIC FOR WENO_S
    const Basis& B;
    
    /// The very Coeffs for one proc
    std::vector<numvector<double, dimS>> SOL;

    /// Full pack of solution coeffs
    std::vector<numvector<double, dimS>> fullSOL;

    /// Pack of solutions in cell centers to export
    std::vector<numvector<double, dimExp>> solToExport;

    /// vector of grad of conservative variables
    std::vector<numvector<double, dimGradCoeff>> S;

    /// Constructor
    Solution(Basis& bas);

    /// Destructor
    ~Solution() {}

    /// Reconstruct solution vector at the point
    numvector<double, dimPh> reconstruct(int iCell, const Point& point) const;

    /// Reconstruct solution component at the point
    double reconstruct(int iCell, const Point& point, Variables var) const;

    /// Reconstruct solution vector using given coeffs
    numvector<double, dimPh> reconstruct(int iCell, const Point& point, const numvector<double, dimS>& SOL) const;

    /// Reconstruct solution component using given coeffs
    double reconstruct(int iCell, const Point& point, const numvector<double, dimS>& SOL, Variables var) const;




};// end Solution

#endif // SOLUTION_H


