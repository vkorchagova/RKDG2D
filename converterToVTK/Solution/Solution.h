#ifndef SOLUTION_H
#define SOLUTION_H

#include <math.h>
#include <functional>
#include <memory>

#include "numvector.h"
#include "defs.h"
#include "Basis.h"
#include "Params.h"


class Solution
{

public:   
    
    //- The very Coeffs
    std::vector<numvector<double, dimS>> SOL;

    //- Reference to the basis
    const Basis& B;

public:

    //- Constructor
    Solution(Basis& bas);

    //- Destructor
    ~Solution() {}

    //- Reconstruct solution at the point
    numvector<double, dimPh> reconstruct(int iCell, const Point& point) const;
    double reconstruct(int iCell, const Point& point, Variables var) const;

    //- Reconstruct solution using given coeffs
    numvector<double, dimPh> reconstruct(int iCell, const Point& point, const numvector<double, dimS>& SOL) const;
    double reconstruct(int iCell, const Point& point, const numvector<double, dimS>& SOL, Variables var) const;

};// end Solution

#endif // SOLUTION_H


