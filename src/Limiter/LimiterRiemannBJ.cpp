#include "LimiterRiemannBJ.h"
#include <omp.h>


using namespace std;

numvector<double, dimS> LimiterRiemannBJ::limitation(const std::vector<std::shared_ptr<Cell>>& stencil)
{ 
    
    numvector<double, dimS> alphaNew;                       // result of limitation

    // p polynoms
    vector<numvector<double, dimS>> pInv;
    vector<numvector<double, dimS>> pNew;

    pInv.reserve(maxPossibleStencilSize);
    pNew.reserve(maxPossibleStencilSize);

    numvector<double, dimPh> solMean;

    int nCells = stencil.size();
    const shared_ptr<Cell>& cell = stencil[0];

    // get p polynoms for each cell

    pInv.resize(nCells);

    //cout << "before sol mean" << endl;
            
    solMean = solution.reconstruct(cell->number,cell->getCellCenter());

    // get pNew polynoms using each edge normal
    for (const shared_ptr<Edge> e : cell->edges)
    {            
        if (e->neibCells.size() == 2)
        {
                // correct normal direction: outside related to troubled cell
            Point n = (( *e->nodes[0] - cell->getCellCenter()) * e->n > 0.0) ? e->n : Point(-e->n);

            L = physics.getL(solMean, n);
            R = physics.getR(solMean, n);

            physics.computeLR(solMean, n, L, R);

            // get Riemann invariants along this direction
            for (size_t k = 0; k < stencil.size(); ++k)
                pInv[k] = conservativeToRiemann(solution.SOL[stencil[k]->number]);

            // limit Riemann invariants
            numvector <double, dimS> pLim = limitChar(stencil, pInv);

            // numvector <double, dimS> res(0.0);

            // for (int i = 0; i < dimPh; ++i)
            // {
            //     res[i*nShapes] =  pLim[0][i*nShapes];
            
            //     for (size_t k = 0; k < nCells; ++k)
            //     {
            //         //res[i*nShapes]     += w[k][i] * p[k][i*nShapes];
            //         res[i*nShapes + 1] += pLim[k][i*nShapes + 1];
            //         res[i*nShapes + 2] += pLim[k][i*nShapes + 2];
            //     }
            // }
                

                /*if (stenc[0]->number == 3163)
                {
                    //cout << "pInv = " << pInv[0] << endl;
                    cout << "pLim = " << pLim << endl;
                    cout << "convert = " << riemannToConservative(pLim, R) << endl;
                }*/

                // project limited solution to conservative variables
            pNew.push_back( riemannToConservative(pLim) );
        }
    }

    double sumArea = 0.0;

    for (int i = 1; i < nCells; ++i)
        sumArea += stencil[i]->getArea();

       /* if (iCell == 3163)
        {
            cout << "sumArea = " << sumArea << endl; 
            cout << "size pNew = " << pNew.size() << endl;  
        }*/

    // construct limited solution
    numvector<double, dimS> res (0.0);
    for (int i = 0; i < nCells-1; ++i)
    {
        res += pNew[i] * (stencil[i+1]->getArea() / sumArea);
            //if (iCell == 3163)
            //    cout << "resInner = " << res << endl; 

    }

          /* if (stenc[0]->number == 18336)
                {
                    cout << "p[0] = " << p[0] << endl;
                    cout << "res = " << res << endl;
                }*/



    alphaNew = res;

    return alphaNew;
}



numvector<double, dimS> LimiterRiemannBJ::conservativeToRiemann
(
    const numvector<double, dimS>& alpha
) const
{
    numvector<double, dimS> res (0.0);

    for (int iSol = 0; iSol < dimPh; ++iSol )
        for (int jSol = 0; jSol < dimPh; ++jSol )
            for (int j = 0; j < nShapes; ++j )
                res[iSol * nShapes + j] +=  L[iSol][jSol] * alpha[jSol * nShapes + j];

    return res;
}

numvector<double, dimS> LimiterRiemannBJ::riemannToConservative
(
    const numvector<double, dimS>& alpha
) const
{
    numvector<double, dimS> res (0);

    for (int iSol = 0; iSol < dimPh; ++iSol )
        for (int jSol = 0; jSol < dimPh; ++jSol )
            for (int j = 0; j < nShapes; ++j )
                res[iSol * nShapes + j] += R[iSol][jSol] * alpha[jSol * nShapes + j];

    return res;

}

numvector<double, dimS> LimiterRiemannBJ::limitChar(const std::vector<std::shared_ptr<Cell>>& stencil, const vector<numvector<double, dimS>>& p)
{
    int nCells = stencil.size();

    // load memory
    std::vector<numvector<double, dimPh>> uMean(nCells);    // mean values
    numvector<double, dimS> alphaNew;                       // result of limitation
    
    // get mean values of linear functions
    uMean.resize(stencil.size());
    for (size_t k = 0; k < nCells; ++k)
        uMean[k] = solution.reconstruct(stencil[k]->number, stencil[k]->getCellCenter(), p[k]);

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
    numvector<double, dimPh> a = getYMin(stencil[0], mI, MI, uMean[0], p[0]);

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

    for (int iSol = 1; iSol < dimPh-2; ++iSol)
        aMin = a[iSol] < aMin ? a[iSol] : aMin;

    // update solution
    alphaNew = p[0];

    for (int iSol = 0; iSol < dimPh-1; ++iSol)
    {
        alphaNew[iSol*nShapes + 1] *= aMin;
        alphaNew[iSol*nShapes + 2] *= aMin;
    }

    alphaNew[(dimPh-1)*nShapes + 1] *= a[0];
    alphaNew[(dimPh-1)*nShapes + 2] *= a[0];

    return alphaNew;
}
