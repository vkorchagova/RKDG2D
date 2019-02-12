#include "LimiterBJ.h"
#include <omp.h>


using namespace std;

numvector<double, dimPh> LimiterBJ::getAlphaL(shared_ptr<Cell>& cell, numvector<double, dimPh>& mI, numvector<double, dimPh>& MI, numvector<double, dimPh>& uMean)
{
    numvector<double, dimPh> y = 0.0;
    numvector<double, dimPh> yMin = 1e9;
    
    //cout << "------" << endl;
    //cout << "Cell No = " << cell->number << endl;
    //for (int e = 0; e < cell->nEntities; ++e)
    for (const shared_ptr<Edge> e : cell->edges)
        for (int i = 0; i < e->nGP; ++i)
        {
            numvector<double, dimPh> pU   = solution.reconstruct(cell->number, e->gPoints[i]);
            numvector<double, dimPh> diff = pU - uMean;
            
        //    cout << "pU - uMean = " << diff << endl;
            
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


void LimiterBJ::limit(vector<numvector<double, dimS>>& alpha)
{
    double ts = omp_get_wtime();

    int n = alpha.size(); // here must be only real cells!
    
    vector<numvector<double, dimS>> alphaNew(n);
    
    int numSol = 0;

    for (int i = 0; i < n; ++i)
        alphaNew[i] = alpha[i];  //no vector = vector???
        
    // mean values

    vector<numvector<double, dimPh>> uMean;

    //use limiter for all cells
    vector<int> troubledCells(n);

    for (int i = 0; i < n; ++i)
        troubledCells[i] = i;


	omp_set_num_threads(NumThreads);
//#pragma omp parallel default(none) \
 shared(myRank, n, alpha, alphaNew, troubledCells) \
 private (numSol, uMean)
//#pragma omp for
    //for (size_t i = 0; i < troubledCells.size(); ++i)
    //for (int iCell : troubledCells)
    for (size_t i = 0; i < n; ++i) // may be in MPI we can iterate through vector without index like Python???
    {
        int iCell = troubledCells[i];

        numvector<double, dimPh> mI =  1e9;
        numvector<double, dimPh> MI = -1e9;


        shared_ptr<Cell> cell = cells[iCell]; 

        // construct list of cells: cell + neighbours

        vector<shared_ptr<Cell>> stenc = { cell }; // the stencil of limitation
        stenc.insert(stenc.end(), cell->neibCells.begin(), cell->neibCells.end());

        int nCells = stenc.size(); // nCells in stencil

        // get mean values of linear functions

        uMean.resize(nCells);

        for (size_t k = 0; k < nCells; ++k)
            uMean[k] = solution.reconstruct(stenc[k]->number, stenc[k]->getCellCenter());
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
        
        // compute a-coeff
        numvector<double, dimPh> a = getAlphaL(cell, mI, MI, uMean[0]);

        alphaNew[iCell] = alpha[iCell];
        
        for (int iSol = 0; iSol < dimPh; ++iSol)
        {
            alphaNew[iCell][iSol*nShapes + 1] *= a[iSol];
            alphaNew[iCell][iSol*nShapes + 2] *= a[iSol];
        }
    }

    alpha = alphaNew;

    
    double te = omp_get_wtime();
    
    //cout << "BJ limiter " << te - ts << endl; 
    lastHope(alpha);
}
