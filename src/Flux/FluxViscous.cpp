#include "FluxViscous.h"
#include <iostream>

using namespace std;

// ------------------ Constructors & Destructors ----------------


// ------------------ Private class methods --------------------

// ------------------ Public class methods ---------------------



numvector<double, dimPh> FluxViscous::evaluate(
        const numvector<double, dimPh>& solLeftRot, const numvector<double, dimPh>& solRightRot,
        const numvector<double, dimGrad>& gradSolLeft, const numvector<double, dimGrad>& gradSolRight
) const
{
    numvector<double, dimPh> fluxOutward = phs.fluxFv(solLeftRot, gradSolLeft);
    numvector<double, dimPh> fluxInward  = phs.fluxFv(solRightRot, gradSolRight);

	//cout << "solLeftRot = " << solLeftRot << "; solRightRot = " << solRightRot << endl;
    //cout << "gradSolLeft = " << gradSolLeft << "; gradSolRight = " << gradSolRight << endl;
    //cout << "fOw = " << fluxOutward << "; fIw = " << fluxInward << endl;

    return 0.5 * (fluxInward + fluxOutward);
}

