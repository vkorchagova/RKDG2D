#include "IndicatorShu.h"

using namespace std;



vector<int> IndicatorShu::checkDiscontinuities() const
{

    vector<int> troubledCells;
    troubledCells.reserve(mesh.nRealCells);

   // cout << "BEG troubledCell.size = " << troubledCells.size() << endl;
    int n = mesh.nRealCells;

    //cout << "Indicator" << endl;

    vector<numvector<double, dimPh>> uMean;
    vector<numvector<double, dimPh>> uMeanOwn;
    uMean.reserve(10);
    uMeanOwn.reserve(10);

    // stencil
    vector<shared_ptr<Cell>> stenc;
    stenc.reserve(10);

//#pragma omp parallel for \
 shared( n, troubledCells) \
 default(none)
//    for (int i = 0; i < n; ++i)
//        troubledCells[i] = i;

#pragma omp parallel for \
            shared(cout, myRank, n, troubledCells) \
            firstprivate(uMean, uMeanOwn, stenc) \
            default(none)
    for (int i = 0; i < n; ++i)
    {
    //check = false;
        int iCell = i;

        const shared_ptr<Cell>& cell = mesh.cells[iCell];

        // construct list of cells: cell + neighbours

        stenc = { cell }; // the stencil of limitation
        stenc.insert(stenc.end(), cell->neibCells.begin(), cell->neibCells.end());

        //stenc = cell->neibCells; // the stencil of limitation
        //stenc.insert(stenc.begin(), cell);


        int nCells = stenc.size(); // nCells in stencil

        // get mean values of linear functions

        uMean.resize(nCells);
        uMeanOwn.resize(nCells);

        for (size_t k = 0; k < nCells; ++k)
        {
           uMean[k] = solution.reconstruct(stenc[k]->number, stenc[0]->getCellCenter());
           uMeanOwn[k] = solution.reconstruct(stenc[k]->number, stenc[k]->getCellCenter());
           //cout << "k = " << k  << endl;
           //cout << stenc[k]->number << endl;
           //cout << stenc[k]->getCellCenter() << endl;
           //cout << "uMean = " << uMean[k] << endl;
           //cout << "uMeanOwn = " << uMeanOwn[k] << endl;

        }

        numvector<double, dimPh> maxFabsPj = -1e9;     // max
        numvector<double, dimPh> sumFabs = 0.0;

        // get maximum from cell averages
        for (size_t k = 0; k < nCells; ++k)
            for (int iSol = 0; iSol < dimPh; ++iSol)
                if (fabs(uMeanOwn[k][iSol]) > maxFabsPj[iSol])
                    maxFabsPj[iSol] = fabs(uMeanOwn[k][iSol]);

        // get smth useful
        for (size_t k = 1; k < nCells; ++k)
            for (int iSol = 0; iSol < dimPh; ++iSol)
                sumFabs[iSol] += fabs(uMean[0][iSol] - uMean[k][iSol]);

        // indicator
        for (size_t k = 0; k < dimPh; ++k)
            if (sumFabs[k] / (maxFabsPj[k] + eps) > Ck)
            {                         
                /*
                cout << "cell#" << iCell << endl;
                for(int i=0;i<nCells;++i)
                    cout << uMean[i] << " ; ";
                cout << endl;
                for(int i=0;i<nCells;++i)
                    cout << uMeanOwn[i] << " ; ";
                cout << endl;
                cout << "sumF: " << sumFabs << endl;
                cout << "maxF: " << maxFabsPj << endl;
                cout << "Ind: ";
                for (size_t k = 0; k < dimPh; ++k)
                    cout <<  (sumFabs[k] / (maxFabsPj[k] + eps)) << ' ';
                cout << endl << endl;
                */
#pragma omp critical
                {
                troubledCells.push_back(iCell);
                //cout << " ind = " << sumFabs[k] / maxFabsPj[k] << ' ';
                //cout << endl;
                }
                break;

            }

    //if (check == true) troubledCells.push_back(iCell);

    }

//#pragma omp parallel for shared (n, troubledCells)

 //  for (size_t i = 0; i < itroubledCells.size(); ++i)
 //     troubledCells.push_back(i);

    return troubledCells;
}
