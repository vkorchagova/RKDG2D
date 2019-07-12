#ifndef FLUX_H
#define FLUX_H

#include "Physics.h"
#include "defs.h"

/// 
/// Abstract class for numerical flux (Riemann solver)
/// 

class Flux
{

public:

    /// Constant reference to physics
    const Physics& phs; // \'-'/

public:

    /// Construct with problem
    Flux(const Physics& phys) : phs(phys) {};

    /// Destructor
    virtual ~Flux() {}

    /// Evaluate numerical flux through one point
    virtual numvector<double, dimPh> evaluate(const numvector<double, dimPh>& solInner, const numvector<double, dimPh>& solOuter) const = 0;

	///
	//// Averaging computations
	///

	/// Compute semisum-averaged sound speed on edge
	double c_av(const numvector<double, dimPh>& solOne, const numvector<double, dimPh>& solTwo) const;

	/// Eigenvalues for X direction ( Roe )
	numvector<double, dimPh> lambdaF_Roe(const numvector<double, dimPh>& solOne, const numvector<double, dimPh>& solTwo) const;

	/// Eigenvalues for X direction ( Einfeldt )
	numvector<double, dimPh> lambdaF_Einfeldt(const numvector<double, dimPh>& solOne, const numvector<double, dimPh>& solTwo) const;

	/// Eigenvalues for X direction ( Toro, pressure-based )
	numvector<double, dimPh> lambdaF_Toro(const numvector<double, dimPh>& solOne, const numvector<double, dimPh>& solTwo) const;

	/// Eigenvalues for X direction ( semisum )
	numvector<double, dimPh> lambdaF_semisum(const numvector<double, dimPh>& solOne, const numvector<double, dimPh>& solTwo) const;

	/// Eigenvalues for X direction
	numvector<double, dimPh> lambdaF(const numvector<double, dimPh>& solOne, const numvector<double, dimPh>& solTwo) const;

};// end Flux

#endif // FLUX_H

