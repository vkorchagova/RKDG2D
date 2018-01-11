#include "LimiterMUSCL.h"

using namespace std;

// --------------------

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

double m(const numvector<double, 3>& slope)
{
    int sign = sgn(slope[0]);
    double minmod = fabs(slope[0]);

    for (int i = 1; i < 3; ++i)
    {
        minmod = (fabs(slope[i]) < minmod) ? fabs(slope[i]) : minmod;

        if (slope[0] * slope[i] < 0.0)
            return 0.0;
    }

    return sign * minmod;
}

// --------------------

void LimiterMUSCL::limit(vector<numvector<double, 5 * nShapes>>& alpha)
{
    //check discontinuities

    problem.setAlpha(alpha);
    vector<int> troubledCells = indicator.checkDiscontinuities();

    // mean values

    numvector<numvector<double, 5>, 3> uMean;

    // limit solution in troubled cells

    for (int iCell : troubledCells)
    {
        // limit in x direction
        // --------------------

        //const Cell& cell = *(indicator.mesh.cells[icell]);

        //vector<shared_ptr<Cell>> neibCells = cell.findNeighbourCellsX();

        numvector<shared_ptr<Cell>, 3> cellsHor = { indicator.mesh.cells[iCell-1], indicator.mesh.cells[iCell], indicator.mesh.cells[iCell+1] };

        // get mean values

        for (int i = 0; i < 3; ++i)
            uMean[i] = cellsHor[i]->reconstructSolution(cellsHor[1]->getCellCenter());

        // limit

        for (int i = 0; i < 5; ++i)
        {
            numvector<double, 3> slope;

            slope[0] = alpha[iCell][i*nShapes + 1] * cellsHor[1]->offsetPhi[1];
            slope[1] = (uMean[i][2] - uMean[i][1]) / cellsHor[1]->h().x() ;
            slope[2] = (uMean[i][1] - uMean[i][0]) / cellsHor[1]->h().x() ;

            alpha[iCell][i*nShapes + 1] = m(slope) / cellsHor[1]->offsetPhi[1];
        }

        problem.setAlpha(alpha);
    }
}

