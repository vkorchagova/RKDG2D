#ifndef PHYSICS_H
#define PHYSICS_H

#include <math.h>
#include <functional>
#include "numvector.h"
#include "Point.h"

#include "defs.h"
#include "Params.h"

///
/// Description of mathematical model: EoS, how to compute natural fluxes, eigvals/eigvectors, speed of sound
///

class Physics
{
     

public:

    /// Heat capacity ratio
    double cpcv;

    /// Covolume constant; 0 for ideal gas
    double covolume;

    /// constant ???
    double lam;

    double Pr;

    double mu;

    /// Constructor
	Physics();

    /// Destructor
    ~Physics();

    /// Calculate pressure using conservative variables
    double getPressure(const numvector<double, dimPh>& sol) const;

    /// Compute sound speed inside cell
    double c(const numvector<double, dimPh>& sol) const;

    /// Compute velocity magnitude
    double magU(const numvector<double, dimPh>& sol) const;

    /// Compute total energy using primitive variables
    double e(double rho, double u, double v, double w, double p) const;

    /// Calculate fluxes in x direction
    numvector<double, dimPh> fluxF(const numvector<double, dimPh>& sol) const;

    /// Calculate fluxes in y direction
    numvector<double, dimPh> fluxG(const numvector<double, dimPh>& sol) const;

    /// Calculate viscous fluxes in x direction
    numvector<double, dimPh> fluxFv(const numvector<double, dimPh>& sol, const numvector<double, dimGrad>& gradSol) const;

    /// Calculate viscous fluxes in y direction
    numvector<double, dimPh> fluxGv(const numvector<double, dimPh>& sol, const numvector<double, dimGrad>& gradSol) const;

    /// Left  eigenvectors
    numvector<numvector<double, dimPh>, dimPh> getL(const numvector<double, dimPh>& sol, const Point& n) const;

	/// Right eigenvectors
    numvector<numvector<double, dimPh>, dimPh> getR(const numvector<double, dimPh>& sol, const Point& n) const;   

};// end Physics

#endif // PHYSICS_H


