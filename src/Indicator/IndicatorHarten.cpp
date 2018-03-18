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

    //check discontinuities on all cells

    for (int i = 0; i < mesh.nCells; ++i)
    {
        // find neighbours

        shared_ptr<Cell> cell = mesh.cells[i];

        vector<shared_ptr<Cell>> neibCellsX = cell->findNeighbourCellsX();
        vector<shared_ptr<Cell>> neibCellsY = cell->findNeighbourCellsY();

        if (neibCellsX.size() == 2)
        {
            // get mean values on the problem cell for x neighbours

            uMean[0] = cell->reconstructSolution(cell->getCellCenter());

            for (size_t k = 1; k < 3; ++k)
                uMean[k] = neibCellsX[k-1]->reconstructSolution(cell->getCellCenter());

            // get coeffs (deriv*coeff near form function)


            // get conditions for trouble cell


            diffSignesForAve = ( (uMean[1] - uMean[0])*(uMean[2] - uMean[0]) < -1e-12 );
            largeDerivRatio = ( fabs(problem.alpha[neibCellsX[0]->number][5]) > kappa*(problem.alpha[i][5]) || \
                                     fabs(problem.alpha[neibCellsX[1]->number][5]) > kappa*(problem.alpha[i][5]) || \
                                     kappa*fabs(problem.alpha[neibCellsX[0]->number][5]) > (problem.alpha[i][5]) || \
                                     kappa*fabs(problem.alpha[neibCellsX[1]->number][5]) > (problem.alpha[i][5]) \
                                   );
        }


        if (diffSignesForAve || largeDerivRatio)
            troubledCells.push_back(i);
    }

//    cout << "\ntroubled cells: " ;
//    for (int iCell : troubledCells)
//        cout << iCell << ' ';

//    cout << endl;

    return troubledCells;
}
