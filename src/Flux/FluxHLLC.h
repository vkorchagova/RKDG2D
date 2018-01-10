#ifndef FLUXHLLC_H
#define FLUXHLLC_H

#include "Flux.h"


class FluxHLLC : public Flux
{
    //- Calculate Star Region velocity
    numvector<double, 5> getUStar (const numvector<double, 5>& sol, double pK, double SK, double cK, double SStar) const;

public:

    //- Construct with problem
    FluxHLLC(const Problem& prb) : Flux(prb) {}

    //- Destructor
    ~FluxHLLC() {}

    //- Evaluate numerical flux through one point
    virtual numvector<double, 5> evaluate(const numvector<double, 5>& solLeft, const numvector<double, 5>& solRight, const Point& n) const override;

};

#endif // FLUXHLLC_H
