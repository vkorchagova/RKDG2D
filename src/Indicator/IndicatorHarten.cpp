#include "IndicatorHarten.h"

using namespace std;

vector<int> IndicatorHarten::checkDiscontinuities() const
{
    vector<int> troubledCells;

    bool diffSignesForAve;
    bool largeDerivRatio;

    // mean values

    vector<numvector<double, 5>> uMean (3);

    // calibration coefficient

    double kappa = 2.0;

    //check discontinuities for all cells

    for (int i = 0; i < mesh.nCells; ++i)
    {
        // find neighbours

        shared_ptr<Cell> cell = mesh.cells[i];

        vector<shared_ptr<Cell>> neibCellsX = cell->findNeighbourCellsX();
        //vector<shared_ptr<Cell>> neibCellsY = cell->findNeighbourCellsY();

        if (neibCellsX.size() == 2)
        {
            // get mean values on the problem cell for x neighbours

            uMean[0] = cell->reconstructSolution(cell->getCellCenter());

            for (size_t k = 1; k < 3; ++k)
                uMean[k] = neibCellsX[k-1]->reconstructSolution(cell->getCellCenter());


            // get conditions for troubled cell


            diffSignesForAve = ( (uMean[1] - uMean[0])*(uMean[2] - uMean[0]) < -1e-7 );

            bool a2 = (fabs(problem.alpha[neibCellsX[0]->number][1]) > kappa*(problem.alpha[i][1]) || \
                    kappa*fabs(problem.alpha[neibCellsX[0]->number][1]) < (problem.alpha[i][1]) );

            bool a3 = (fabs(problem.alpha[neibCellsX[1]->number][1]) > kappa*(problem.alpha[i][1]) || \
                    kappa*fabs(problem.alpha[neibCellsX[1]->number][1]) < (problem.alpha[i][1]) );

            largeDerivRatio = ( a2 || a3);
        }


        if (diffSignesForAve && largeDerivRatio)
            troubledCells.push_back(i);
    }

    cout << "\ntroubled cells: " ;
    for (int iCell : troubledCells)
        cout << iCell << ' ';

    cout << endl;

    return troubledCells;
}
