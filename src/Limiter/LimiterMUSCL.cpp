#include "LimiterMUSCL.h"

using namespace std;

// --------------------

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

double m(const numvector<double, 3>& slope)
{
    double sign = sgn(slope[0]);
    double minmod = fabs(slope[0]);

    for (int i = 1; i < 3; ++i)
    {
        minmod = (fabs(slope[i]) < minmod) ? fabs(slope[i]) : minmod;

        if (slope[0] * slope[i] < 0)
            return 0.0;
    }

    return sign * minmod;
}

// --------------------

void LimiterMUSCL::limit(vector<numvector<double, 5 * nShapes>>& alpha)
{
    //check discontinuities

    problem.setAlpha(alpha);
    vector<double> ind = indicator.checkDiscontinuities();



    // initialize needed arrays

    numvector<numvector<double, 5>, 3> uMean;

    // limit solution in troubled cells

    for (size_t icell = 0; icell < alpha.size(); ++icell)
    {
        if (ind[icell] > 1.0)
        {

            cout << "troubled cell #" << icell << " " << ind[icell] << endl;

            // limit in x direction
            // --------------------

            //const Cell& cell = *(indicator.mesh.cells[icell]);

            //vector<shared_ptr<Cell>> neibCells = cell.findNeighbourCellsX();

            numvector<shared_ptr<Cell>, 3> cellsHor = { indicator.mesh.cells[icell-1], indicator.mesh.cells[icell], indicator.mesh.cells[icell+1] };

            // get mean values

            for (int i = 0; i < 3; ++i)
                uMean[i] = cellsHor[i]->reconstructSolution(cellsHor[1]->getCellCenter());

            // limit

            for (int i = 0; i < 5; ++i)
            {
                numvector<double, 3> slope;

                slope[0] = alpha[icell][i*nShapes + 1];
                slope[1] = (uMean[i][2] / cellsHor[2]->offsetPhi[1]  - uMean[i][1] / cellsHor[1]->offsetPhi[1]) / 0.5 / (cellsHor[2]->h().x() + cellsHor[1]->h().x());
                slope[2] = (uMean[i][1] / cellsHor[1]->offsetPhi[1]  - uMean[i][0] / cellsHor[0]->offsetPhi[1]) / 0.5 / (cellsHor[1]->h().x() + cellsHor[0]->h().x());

                alpha[icell][i*nShapes + 1] = m(slope);
            }

        }
    }

    return;
}

