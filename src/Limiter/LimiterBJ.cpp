#include "LimiterBJ.h"
#include <omp.h>


using namespace std;

LimiterBJ::LimiterBJ(
        const Mesh& mesh, 
        Solution& sln,
        const Physics& phs,
        const Indicator& ind) : Limiter(mesh, sln, phs, ind) 
{
    //uMean.reserve(maxPossibleStencilSize);
}

numvector<double, dimPh> LimiterBJ::getYMin(
    const shared_ptr<Cell>& cell, 
    const numvector<double, dimPh>& mI, 
    const numvector<double, dimPh>& MI, 
    const numvector<double, dimPh>& uMean)
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


vector<shared_ptr<Cell>> LimiterBJ::getStencilFor(const std::shared_ptr<Cell>& cell)
{
    vector<shared_ptr<Cell>> stencil = { cell }; 
    stencil.insert(stencil.end(), cell->neibCells.begin(), cell->neibCells.end());

    return stencil;
}


numvector<double, dimS> LimiterBJ::limitation(const std::vector<std::shared_ptr<Cell>>& stencil)
{
    int nCells = stencil.size();

    // load memory
    std::vector<numvector<double, dimPh>> uMean(nCells);    // mean values
    numvector<double, dimS> alphaNew;                       // result of limitation
    
    // get mean values of linear functions
    uMean.resize(stencil.size());
    for (size_t k = 0; k < nCells; ++k)
        uMean[k] = solution.reconstruct(stencil[k]->number, stencil[k]->getCellCenter());

    // get minimum from cell averages
    numvector<double, dimPh> mI = 1e9;      // min
    numvector<double, dimPh> MI = -1e9;     // max
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
    numvector<double, dimPh> a = getYMin(stencil[0], mI, MI, uMean[0]);

    // update solution
    alphaNew = solution.SOL[stencil[0]->number];

    for (int iSol = 0; iSol < dimPh; ++iSol)
    {
        alphaNew[iSol*nShapes + 1] *= a[iSol];
        alphaNew[iSol*nShapes + 2] *= a[iSol];
    }

    return alphaNew;
}
