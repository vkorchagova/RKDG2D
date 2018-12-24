#include "BoundarySlip.h"
#include "Basis.h"

void BoundarySlip::applyBoundary(Basis& basis) const
{
    for (int i = 0; i < patch.edgeGroup.size(); ++i)
    {
        int numReal    = patch.edgeGroup[i]->neibCells[0]->number;
        int numFiction = patch.edgeGroup[i]->neibCells[1]->number;

        basis.phiCoeffs[numFiction] = basis.phiCoeffs[numReal];

        basis.phiCoeffs[numFiction][1] = 0.0;
        basis.phiCoeffs[numFiction][2] = 0.0;
    }
}
