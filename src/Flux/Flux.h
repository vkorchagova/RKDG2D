#ifndef FLUX_H
#define FLUX_H

#include "Physics.h"
#include "defs.h"

/// 
/// Abstract class for numerical flux (Riemann solver)
/// 

class Flux
{
private:



protected:

    /// Constant reference to physics
    const Physics& phs;

    ///
    //// Averaging computations
    ///

    /// Compute semisum-averaged sound speed on edge
    double c_av(const numvector<double, dimPh>& solOne, const numvector<double, dimPh>& solTwo) const;
 
    /// Eigenvalues for X direction ( Roe ) // UPDATE FOR COVOLUME
    numvector<double, dimPh> lambdaF_Roe(const numvector<double, dimPh>& solOne, const numvector<double, dimPh>& solTwo) const;

    /// Eigenvalues for X direction ( Einfeldt ) // UPDATE FOR COVOLUME
    numvector<double, dimPh> lambdaF_Einfeldt(const numvector<double, dimPh>& solOne, const numvector<double, dimPh>& solTwo) const;

    /// Eigenvalues for X direction ( Toro, pressure-based ) // UPDATE FOR COVOLUME
    numvector<double, dimPh> lambdaF_Toro(const numvector<double, dimPh>& solOne, const numvector<double, dimPh>& solTwo) const;

    /// Eigenvalues for X direction ( semisum ) // UPDATE FOR COVOLUME
    numvector<double, dimPh> lambdaF_semisum(const numvector<double, dimPh>& solOne, const numvector<double, dimPh>& solTwo) const;

    /// Eigenvalues for X direction // UPDATE FOR COVOLUME
    numvector<double, dimPh> lambdaF(const numvector<double, dimPh>& solOne, const numvector<double, dimPh>& solTwo) const;

public:

    /// Constructor
    //Flux(const Physics& phys) : phs(phys) {};
    Flux(const Physics& phys);

    /// Destructor
    virtual ~Flux() {}

    /// Evaluate numerical flux through one point
    //virtual numvector<double, dimPh> evaluate(const numvector<double, dimPh>& solInner, const numvector<double, dimPh>& solOuter) const = 0;

    virtual numvector<double, dimPh> evaluate(const numvector<double, dimPh>& solInner, const numvector<double, dimPh>& solOuter,
                                              const numvector<double, dimGrad>& gradSolInner = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                              const numvector<double, dimGrad>& gradSolOuter = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}) const = 0;

};// end Flux

#endif // FLUX_H

