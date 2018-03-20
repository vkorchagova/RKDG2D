#include "IndicatorHarten.h"

using namespace std;

vector<int> IndicatorHarten::checkDiscontinuities() const
{
    vector<int> troubledCells;

    bool diffSignesForAve;
    bool largeDerivRatio;

    // mean values

    vector<numvector<double, 5>> uMean (3);
    double meanCondition;
    // calibration coefficient

    double kappa = 2.0;

    bool a2, a3;

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

            meanCondition = (uMean[1][0] - uMean[0][0])*(uMean[2][0] - uMean[0][0]);
            diffSignesForAve = ( meanCondition < -1e-7 );

            // neibCellsX[0]->number = number of left neighbour
            // neibCellsX[1]->number = number of right neighbour
            // i = number of considered cell
            // [1] = coefficient for density for 1st basis function

            a2 = (fabs(problem.alpha[neibCellsX[0]->number][1]) > kappa*fabs(problem.alpha[i][1]) || \
                    kappa*fabs(problem.alpha[neibCellsX[0]->number][1]) < fabs(problem.alpha[i][1]) );

            a3 = (fabs(problem.alpha[neibCellsX[1]->number][1]) > kappa*fabs(problem.alpha[i][1]) || \
                    kappa*fabs(problem.alpha[neibCellsX[1]->number][1]) < fabs(problem.alpha[i][1]) );

            largeDerivRatio = ( a2 || a3);
        }


        if (diffSignesForAve && largeDerivRatio)
        {
//            cout << "cell number = " << i << endl;
//            cout << "mean condition: " << meanCondition << endl;
//            for (int jj = 0; jj < 3; ++jj)
//                cout << uMean[jj][0] << endl;
//            //cout << "1st condition:" << (uMean[1][0] - uMean[0][0])*(uMean[2][0] - uMean[0][0]) << endl;
//            cout << "v: "<< fabs(problem.alpha[i][1]) << endl;
//            cout << "v3: "<< fabs(problem.alpha[neibCellsX[0]->number][1]) << endl;
//            cout << "v4: "<< fabs(problem.alpha[neibCellsX[1]->number][1]) << endl;

//            cout << "a1 = " << diffSignesForAve << ", a2 = "<< a2 << " a3 = " << a3 << endl;

            troubledCells.push_back(i);
        }
    }

    cout << "-----\ntroubled cells: " ;
    for (int iCell : troubledCells)
        cout << iCell << ' ';

    cout << endl;
     cout << "-----\n";

    return troubledCells;
}
