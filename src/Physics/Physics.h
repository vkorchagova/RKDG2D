#ifndef PHYSICS_H
#define PHYSICS_H

#include <math.h>
#include <functional>
#include "numvector.h"
#include "Point.h"

#include "defs.h"
#include "Params.h"

enum Variables { RHO, RHOU, RHOV, RHOW, E};

class Physics
{

public:

    //- Heat capacity ratio
    double cpcv;

    //- Constructor
	Physics();
    //- Destructor
    ~Physics();

    //- Calculate pressure using conservative variables
    double getPressure(const numvector<double, dimPh>& sol) const;

    //- Compute sound speed inside cell
    double c(const numvector<double, dimPh>& sol) const;

    //- Calculate fluxes in x direction
    numvector<double, dimPh> fluxF(const numvector<double, dimPh>& sol) const;

    //- Calculate fluxes in y direction
    numvector<double, dimPh> fluxG(const numvector<double, dimPh>& sol) const;

    //- Left  eigenvectors
    //std::pair<numvector<numvector<double, dimPh>, dimPh>, numvector<numvector<double, dimPh>, dimPh>> getL(const numvector<double, dimPh>& sol) const;
	numvector<numvector<double, dimPh>, dimPh> getL(const numvector<double, dimPh>& sol, const Point& n) const;
	numvector<numvector<double, dimPh>, dimPh> getLx(const numvector<double, dimPh>& sol) const;
	numvector<numvector<double, dimPh>, dimPh> getLy(const numvector<double, dimPh>& sol) const;
	//- Right eigenvectors
    //std::pair<numvector<numvector<double, dimPh>, dimPh>, numvector<numvector<double, dimPh>, dimPh>> getR(const numvector<double, dimPh>& sol) const;
	numvector<numvector<double, dimPh>, dimPh> getR(const numvector<double, dimPh>& sol, const Point& n) const;
	numvector<numvector<double, dimPh>, dimPh> getRx(const numvector<double, dimPh>& sol) const;
	numvector<numvector<double, dimPh>, dimPh> getRy(const numvector<double, dimPh>& sol) const;

};// end Physics

#endif // PHYSICS_H


