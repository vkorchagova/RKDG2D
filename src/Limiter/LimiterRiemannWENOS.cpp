#include "LimiterRiemannWENOS.h"

using namespace std;

numvector<double, dimS> LimiterRiemannWENOS::limitation(const std::vector<std::shared_ptr<Cell>>& stencil)
{ 
    
    numvector<double, dimS> alphaNew;                       // result of limitation

    // p polynoms
    vector<numvector<double, dimS>> pInv;
    vector<numvector<double, dimS>> pNew;

    pInv.reserve(maxPossibleStencilSize);
    pNew.reserve(maxPossibleStencilSize);


    numvector<numvector<double, dimPh>, dimPh>  L;
    numvector<numvector<double, dimPh>, dimPh>  R;
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

            // get Riemann invariants along this direction
            for (size_t k = 0; k < stencil.size(); ++k)
                pInv[k] = conservativeToRiemann(solution.SOL[stencil[k]->number], L);

            // limit Riemann invariants
            numvector <double, dimS> pLim = limitP(stencil, pInv);

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
            pNew.push_back( riemannToConservative(pLim, R) );
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



numvector<double, dimS> LimiterRiemannWENOS::conservativeToRiemann
(
    const numvector<double, dimS>& alpha, 
    const numvector<numvector<double, dimPh>, dimPh>& L
) const
{
    numvector<double, dimS> res (0.0);

    for (int iSol = 0; iSol < dimPh; ++iSol )
        for (int jSol = 0; jSol < dimPh; ++jSol )
            for (int j = 0; j < nShapes; ++j )
                res[iSol * nShapes + j] +=  L[iSol][jSol] * alpha[jSol * nShapes + j];

    return res;
}

numvector<double, dimS> LimiterRiemannWENOS::riemannToConservative
(
    const numvector<double, dimS>& alpha, 
    const numvector<numvector<double, dimPh>, dimPh>& R
) const
{
    numvector<double, dimS> res (0);

    for (int iSol = 0; iSol < dimPh; ++iSol )
        for (int jSol = 0; jSol < dimPh; ++jSol )
            for (int j = 0; j < nShapes; ++j )
                res[iSol * nShapes + j] += R[iSol][jSol] * alpha[jSol * nShapes + j];

    return res;

}

vector<shared_ptr<Cell>> LimiterRiemannWENOS::getStencilFor(const std::shared_ptr<Cell>& cell)
{
    vector<shared_ptr<Cell>> stencil = { cell }; 
    stencil.insert(stencil.end(), cell->neibCells.begin(), cell->neibCells.end());
    return stencil;
}



