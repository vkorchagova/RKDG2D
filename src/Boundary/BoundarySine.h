#ifndef BOUNDARYSINE_H
#define BOUNDARYSINE_H

#include "numvector.h"
#include "TimeClass.h"
#include "Boundary.h"

class Problem;


class BoundarySine : public Boundary
{

protected:

    //- Amplitude
    double a_;

    //- Frequency
    double f_;

    //- Reference value
    numvector<double,5> u0_;

    //- Reference to time
    const Time& time_;

    //- Reference to problem
    const Problem& problem_;

public:

    //- Constructor
    BoundarySine(double a, double f, const Time& t, const Problem& prb, const numvector<double,5>& u0 = {1,0,0,0,1}) : Boundary(), a_(a), f_(f), u0_(u0), time_(t), problem_(prb) { type = "sine";}

    //- Set ref value
    void setRefValue(const numvector<double,5>& u0) {u0_ = u0;}

    //- Apply boundary condition
    virtual numvector<double, 5> applyBoundary(const numvector<double, 5>& solLeft = {0.0, 0.0, 0.0, 0.0, 0.0}, const Point& n = Point({0.0,0.0})) const override;

};

#endif // BOUNDARYSINE_H
