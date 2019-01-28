#include "BoundarySlip.h"
#include "Basis.h"

void BoundarySlip::applyBoundary(std::vector<numvector<double, dimS>>& coeffs) const
{
    for (int i = 0; i < patch.edgeGroup.size(); ++i)
    {
        int numReal    = patch.edgeGroup[i]->neibCells[0]->number;
        int numFiction = patch.edgeGroup[i]->neibCells[1]->number;

        coeffs[numFiction] = coeffs[numReal];

        // for (int iSol = 0; iSol < dimPh; ++iSol)
        // {
        //     coeffs[numFiction][iSol * nShapes + 1] = 0.0;
        //     coeffs[numFiction][iSol * nShapes + 2] = 0.0;
        // }
        
    }
}

numvector<double, dimPh> BoundarySlip::getSolOuter (const numvector<double, 5>& solLeft) const
{
    return { solLeft[0], -solLeft[1], solLeft[2], solLeft[3], solLeft[4]};
}