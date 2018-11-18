#ifndef BOUNDARYSINEDIR_H
#define BOUNDARYSINEDIR_H

#include "BoundarySine.h"

class BoundarySineDir : public BoundarySine
{
public:

    BoundarySineDir(double a, double f, const Time& t, const Problem& prb, const numvector<double,5>& u0 = {1,0,0,0,1}) : BoundarySine(a, f, t, prb, u0) { type = "sineDir";}

    //- Apply boundary condition
    virtual numvector<double, 5> applyBoundary(const numvector<double, 5>& solLeft = {0.0, 0.0, 0.0, 0.0, 0.0}, const Point& n = Point({1.0,0.0})) const override;

};

#endif // BOUNDARYSINEDIR_H
