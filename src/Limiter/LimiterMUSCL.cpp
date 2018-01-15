#include "LimiterMUSCL.h"

using namespace std;

// --------------------

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

double m(const vector<double>& slope)
{
    int sign = sgn(slope[0]);
    double minmod = fabs(slope[0]);


    for (size_t i = 1; i < slope.size(); ++i)
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

    vector<numvector<double, 5>> uMeanX;
    vector<numvector<double, 5>> uMeanY;

    // limit solution in troubled cells

    for (int iCell : troubledCells)
    {
        shared_ptr<Cell> cell = indicator.mesh.cells[iCell];

        vector<shared_ptr<Cell>> neibCellsX = cell->findNeighbourCellsX();
        vector<shared_ptr<Cell>> neibCellsY = cell->findNeighbourCellsY();

        vector<shared_ptr<Cell>> cellsHor = { cell };
        vector<shared_ptr<Cell>> cellsVer = { cell };

        cellsHor.insert(cellsHor.end(), neibCellsX.begin(), neibCellsX.end());
        cellsVer.insert(cellsVer.end(), neibCellsY.begin(), neibCellsY.end());

        //numvector<shared_ptr<Cell>, 3> cellsHor = { indicator.mesh.cells[iCell-1], indicator.mesh.cells[iCell], indicator.mesh.cells[iCell+1] };

        // get mean values

        for (size_t i = 0; i < cellsHor.size(); ++i)
            uMeanX.push_back( cellsHor[i]->reconstructSolution(cellsHor[0]->getCellCenter()) );

        for (size_t i = 0; i < cellsVer.size(); ++i)
            uMeanY.push_back( cellsVer[i]->reconstructSolution(cellsVer[0]->getCellCenter()) );

        // limit

        for (int i = 0; i < 5; ++i)
        {
            vector<double> slopeX;
            vector<double> slopeY;

            slopeX.push_back(alpha[iCell][i*nShapes + 1]);
            slopeY.push_back(alpha[iCell][i*nShapes + 2]);

            for (size_t j = 1; j < cellsHor.size(); ++j)
                slopeX.push_back( sgn(cellsHor[j]->getCellCenter().x() - cellsHor[0]->getCellCenter().x()) * (uMeanX[i][j] - uMeanX[i][0]) / cellsHor[0]->h().x() / cellsHor[0]->offsetPhi[1] );

            for (size_t j = 1; j < cellsVer.size(); ++j)
                slopeY.push_back( sgn(cellsVer[j]->getCellCenter().y() - cellsVer[0]->getCellCenter().y()) * (uMeanY[i][j] - uMeanY[i][0]) / cellsVer[0]->h().y() / cellsHor[0]->offsetPhi[2] );

            //cout << "i = "<< i << ' ' << m(slopeX) / cellsHor[0]->offsetPhi[1] << ' ' << m(slopeY) / cellsVer[0]->offsetPhi[2] << endl;

            alpha[iCell][i*nShapes + 1] = m(slopeX);
            alpha[iCell][i*nShapes + 2] = m(slopeY);
        }

        problem.setAlpha(alpha);
    }
}

