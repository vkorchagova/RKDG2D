#ifndef BOUNDARYCONSTANT_H
#define BOUNDARYCONSTANT_H

#include "Boundary.h"

class BoundaryConstant : public Boundary
{
    numvector<double, 5> fixedValue;

public:

    //- Constructor
    BoundaryConstant(const numvector<double, 5>& defValue) : fixedValue(defValue) {};

    //- Apply boundary condition
    numvector<double, 5> applyBoundary(const numvector<double, 5>& solLeft = {0.0, 0.0, 0.0, 0.0, 0.0}) const override;

};

#endif // BOUNDARYCONSTANT_H
