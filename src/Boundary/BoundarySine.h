#ifndef BOUNDARYSINE_H
#define BOUNDARYSINE_H

#include "numvector.h"
#include "Time.h"
#include "Boundary.h"


class BoundarySine : public Boundary
{
    //- Amplitude
    double a_;

    //- Frequency
    double f_;

    //- Reference value
    numvector<double,5> u0_;

    //- Reference to time
    const Time& time_;

public:

    //- Constructor
    BoundarySine(double a, double f, const Time& t, const numvector<double,5>& u0 = {0,0,0,0,0}) : a_(a), f_(f), u0_(u0), time_(t) {}

    //- Set ref value
    void setRefValue(const numvector<double,5>& u0) {u0_ = u0;}

    //- Apply boundary condition
    numvector<double, 5> applyBoundary(const numvector<double, 5>& solLeft = {0.0, 0.0, 0.0, 0.0, 0.0}) const override;

};

#endif // BOUNDARYSINE_H
