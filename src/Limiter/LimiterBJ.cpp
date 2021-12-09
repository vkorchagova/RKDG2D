#include "LimiterBJ.h"
#include <omp.h>


using namespace std;

LimiterBJ::LimiterBJ(
        const Mesh& mesh, 
        Solution& sln,
        const Physics& phs,
        const Indicator& ind,
        Buffers& _buf) : Limiter(mesh, sln, phs, ind, _buf) 
{
    //uMean.reserve(maxPossibleStencilSize);
}

numvector<double, dimPh> LimiterBJ::getYMin(
    const shared_ptr<Cell>& cell, 
    const numvector<double, dimPh>& mI, 
    const numvector<double, dimPh>& MI, 
    const numvector<double, dimPh>& uMean,
    const numvector<double, dimS>& p)
{
    numvector<double, dimPh> y = 0.0;
    numvector<double, dimPh> yMin = 1e9;
    
    //cout << "------" << endl;
    //cout << "Cell No = " << cell->number << endl;
    //for (int e = 0; e < cell->nEntities; ++e)
    for (const shared_ptr<Edge> e : cell->edges)
        for (int i = 0; i < e->nGP; ++i)
        {
            numvector<double, dimPh> pU   = solution.reconstruct(cell->number, e->gPoints[i], p);
            numvector<double, dimPh> diff = pU - uMean;
            
        //    cout << "pU - uMean = " << diff << endl;

            // numvector<double, dimPh> physQuite = {1.204, 1.0, 1.0, 1.0, 253312};
            
            for (int i = 0; i < dimPh; ++i)
            {
                if (diff[i] > 1e-3 * max(1.0, fabs(uMean[i])))
                    y[i] = (MI[i] - uMean[i]) / diff[i];
                else if (diff[i] < -1e-3 * max(1.0, fabs(uMean[i])))
                    y[i] = (mI[i] - uMean[i]) / diff[i];
                else
                    y[i] = 1.0;
            
                // yMin[i] = y[i] < yMin[i] ? y[i] : yMin[i];
                // yMin[i] = yMin[i] > 1.0 ? 1.0 : yMin[i];

                // michalak
                double yStar = 1.5;
                double yRel = y[i] / yStar;
                double yCur = y[i] < 1.0 ? y[i] + (3.0 - 2.0 * yStar) * yRel * yRel + (yStar - 2.0) * yRel * yRel * yRel : 1.0;
                yMin[i] = yCur < yMin[i] ? yCur : yMin[i];
                
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
    numvector<double, dimPh> a = getYMin(stencil[0], mI, MI, uMean[0], solution.SOL[stencil[0]->number]);

    // // double aMin = 1e+9;

    // // for (int iSol = 0; iSol < dimPh; ++iSol)
    // //     aMin = a[iSol] < aMin ? a[iSol] : aMin;

    // // update solution
    // alphaNew = solution.SOL[stencil[0]->number];

    // alphaNew[1] *= a[0];
    // alphaNew[2] *= a[0];

    // for (int iSol = 1; iSol < dimPh; ++iSol)
    // {
    //     alphaNew[iSol*nShapes + 1] = alphaNew[1] * uMean[0][iSol] / uMean[0][0]; //a[0];
    //     alphaNew[iSol*nShapes + 2] = alphaNew[2] * uMean[0][iSol] / uMean[0][0]; //a[0];
    // }

    double aMin = 1e+9;

    for (int iSol = 1; iSol < 3; ++iSol)
        aMin = a[iSol] < aMin ? a[iSol] : aMin;

    // update solution
    alphaNew = solution.SOL[stencil[0]->number];

    for (int iSol = 1; iSol < 3; ++iSol)
    {
        alphaNew[iSol*nShapes + 1] *= a[0];
        alphaNew[iSol*nShapes + 2] *= a[0];
    }
    
    alphaNew[1] *= a[0];
    alphaNew[2] *= a[0];

    alphaNew[(dimPh-1)*nShapes + 1] *= a[0];
    alphaNew[(dimPh-1)*nShapes + 2] *= a[0];

    return alphaNew;
}
