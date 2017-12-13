#ifndef FLUX_H
#define FLUX_H

#include "Problem.h"

class Flux
{

public:
    Problem *problem;

public:

    //- Default constructor
    Flux() {}

    //- Construct with problem
    Flux(Problem &prb);

    //- Destructor
    virtual ~Flux() {}
    
    //- Overloaded "=" operator
    Flux& operator=(const Flux& flx) { problem = flx.problem; return *this;}

    //- Evaluate numerical flux through one point
    virtual numvector<double, 5> evaluate(const numvector<double, 5>& solL, const numvector<double, 5>& solR) = 0;

};// end Flux

#endif // FLUX_H

