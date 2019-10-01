#include "IndicatorBJ.h"

using namespace std;


numvector<double, dimPh> IndicatorBJ::getYMin(const shared_ptr<Cell>& cell, const numvector<double, dimPh>& mI, const numvector<double, dimPh>& MI, const  numvector<double, dimPh>& uMean) const
{
    numvector<double, dimPh> y = 0.0;
    numvector<double, dimPh> yMin = 1e9;

    for (const shared_ptr<Edge> e : cell->edges)
        for (int i = 0; i < e->nGP; ++i)
        {
            numvector<double, dimPh> pU   = solution.reconstruct(cell->number, e->gPoints[i]);
            numvector<double, dimPh> diff = pU - uMean;

            for (int i = 0; i < dimPh; ++i)
            {
                if (diff[i] > 1e-6)
                    y[i] = (MI[i] - uMean[i]) / diff[i];
                else if (diff[i] < -1e-6)
                    y[i] = (mI[i] - uMean[i]) / diff[i];
                else
                    y[i] = 1.0;

                yMin[i] = y[i] < yMin[i] ? y[i] : yMin[i];
                yMin[i] = yMin[i] > 1.0 ? 1.0 : yMin[i];

                if (yMin[i] < 0)
                    cout << "yMin < 0 --- so strange!" << endl;
            }
        }

    //cout << "yMin = " << yMin << endl;

    return yMin;
}


vector<int> IndicatorBJ::checkDiscontinuities() const
{

    vector<int> troubledCells;
    troubledCells.reserve(mesh.nRealCells);

   // cout << "BEG troubledCell.size = " << troubledCells.size() << endl;
    int n = mesh.nRealCells;

    //cout << "Indicator" << endl;

    vector<numvector<double, dimPh>> uMean;
    uMean.reserve(10);

    // stencil
    vector<shared_ptr<Cell>> stenc;
    stenc.reserve(10);

//#pragma omp parallel for \
 shared( n, troubledCells) \
 default(none)
//    for (int i = 0; i < n; ++i)
//        troubledCells[i] = i;

#pragma omp parallel for \
            shared(myRank, n, troubledCells) \
            firstprivate(uMean, stenc) \
            default(none)
    for (int i = 0; i < n; ++i)
    {
    //check = false;
        int iCell = i;

        numvector<double, dimPh> mI = 1e9;
        numvector<double, dimPh> MI = -1e9;

        const shared_ptr<Cell>& cell = mesh.cells[iCell];

        // construct list of cells: cell + neighbours

        stenc = { cell }; // the stencil of limitation
        stenc.insert(stenc.end(), cell->neibCells.begin(), cell->neibCells.end());

        //stenc = cell->neibCells; // the stencil of limitation
        //stenc.insert(stenc.begin(), cell);


        int nCells = stenc.size(); // nCells in stencil

        // get mean values of linear functions

        uMean.resize(nCells);

        for (size_t k = 0; k < nCells; ++k)
        {
//    	cout << "k = " << k  << endl;
//    	cout << stenc[k]->number << endl;
//    	cout << stenc[k]->getCellCenter() << endl;
//    	cout << "uMean = " << uMean[k] << endl;
           uMean[k] = solution.reconstruct(stenc[k]->number, stenc[k]->getCellCenter());
        }
        // here tooooo strange to reconstruct solution through number... may be optimal way exists?..


        // get minimum from cell averages
        for (size_t k = 0; k < nCells; ++k)
            for (int iSol = 0; iSol < dimPh; ++iSol)
                if (uMean[k][iSol] < mI[iSol])
                 mI[iSol] = uMean[k][iSol];

        // get maximum from cell averages
        for (size_t k = 0; k < nCells; ++k)
            for (int iSol = 0; iSol < dimPh; ++iSol)
                if (uMean[k][iSol] > MI[iSol])
                    MI[iSol] = uMean[k][iSol];

        // indicator
        numvector<double, dimPh> a = getYMin(cell, mI, MI, uMean[0]);

        for (size_t k = 0; k < dimPh; ++k)
            if (a[k] < 1.0 && a[k]>1e-6)
            {
#pragma omp critical
                troubledCells.push_back(iCell);
                break;
            }

    //if (check == true) troubledCells.push_back(iCell);

    }

//#pragma omp parallel for shared (n, troubledCells)

 //  for (size_t i = 0; i < itroubledCells.size(); ++i)
 //     troubledCells.push_back(i);

    return troubledCells;
}
