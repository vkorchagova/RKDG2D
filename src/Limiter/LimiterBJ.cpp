#include "LimiterBJ.h"
#include <omp.h>

using namespace std;

numvector<double, 5> getAlphaL(shared_ptr<Cell>& cell, numvector<double, 5>& mI, numvector<double, 5>& MI, numvector<double, 5>& uMean)
{
    numvector<double, 5> y = 0.0;
    numvector<double, 5> yMin = 1e9;
    
    //cout << "------" << endl;
    //cout << "Cell No = " << cell->number << endl;
    for (int e = 0; e < cell->nEntities; ++e)
    for (int i = 0; i < cell->edges[e]->nGP; ++i)
    {
        numvector<double, 5> pU   = cell->reconstructSolution(cell->edges[e]->gPoints[i]);
        numvector<double, 5> diff = pU - uMean;
        
    //    cout << "pU - uMean = " << diff << endl;
        
        for (int i = 0; i < 5; ++i)
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


void LimiterBJ::limit(vector<numvector<double, 5 * nShapes>>& alpha)
{
    double ts = omp_get_wtime();
    int n = alpha.size();
    
    problem.setAlpha(alpha);
    
    vector<numvector<double, 5 * nShapes>> alphaNew(n);
    
    int numSol = 0;

#pragma omp parallel for \
shared(n, alphaNew, alpha) \
default(none)
    for (int i = 0; i < n; ++i)
        alphaNew[i] = alpha[i];
        
    // mean values

    vector<numvector<double,5>> uMean;


        // limit solution in troubled cells
#pragma omp parallel for \
shared(alpha, alphaNew,numSol) \
private(uMean) \
default(none)
        //for (size_t i = 0; i < troubledCells.size(); ++i)
        //for (int iCell : troubledCells)
        for (int iCell = 0; iCell < indicator.mesh.nCells; ++iCell)
        {
            numvector<double, 5> mI =  1e9;
            numvector<double, 5> MI = -1e9;

            //int iCell = troubledCells[i];

            shared_ptr<Cell> cell = indicator.mesh.cells[iCell];

            // construct list of cells: cell + neighbours

            vector<shared_ptr<Cell>> cells = { cell };
            cells.insert(cells.end(), cell->neibCells.begin(), cell->neibCells.end());

            int nCells = cells.size();

            // get mean values of linear functions

            uMean.resize(nCells);

            for (size_t k = 0; k < nCells; ++k)
                uMean[k] = cells[k]->reconstructSolution(cells[k]->getCellCenter());

            // get minimum from cell averages
            for (size_t k = 0; k < nCells; ++k)
                for (int iSol = 0; iSol < 5; ++iSol)
                    if (uMean[k][iSol] < mI[iSol])
                        mI[iSol] = uMean[k][iSol];

            // get maximum from cell averages
            for (size_t k = 0; k < nCells; ++k)
                for (int iSol = 0; iSol < 5; ++iSol)
                    if (uMean[k][iSol] > MI[iSol])
                        MI[iSol] = uMean[k][iSol];
            
            // compute a-coeff
            numvector<double, 5> a = getAlphaL(cell, mI, MI, uMean[0]);

            alphaNew[iCell] = alpha[iCell];
            
            for (int iSol = 0; iSol < 5; ++iSol)
            {
                alphaNew[iCell][iSol*nShapes + 1] *= a[iSol];
                alphaNew[iCell][iSol*nShapes + 2] *= a[iSol];
            }
        }

    alpha = alphaNew;

    
    double te = omp_get_wtime();
    
    //cout << "BJ limiter " << te - ts << endl; 
}
