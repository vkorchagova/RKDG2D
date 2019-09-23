#include "LimiterBJVertex.h"
#include <omp.h>


using namespace std;

LimiterBJVertex::LimiterBJVertex(
        const Mesh& mesh, 
        Solution& sln,
        const Physics& phs,
        const Indicator& ind) : Limiter(mesh, sln, phs, ind) 
{
    // uMean.reserve(maxPossibleStencilSize);
}


numvector<double, dimPh> LimiterBJVertex::getYMin(
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
    {
        // check all Gauss points
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

        // check all vertices
        for (auto c : cell->neibCellsVertex )
            for (auto n : c->nodes)
            {
                if (n == e->nodes[1])
                {
                    numvector<double, dimPh> pU   = solution.reconstruct(cell->number, *(e->nodes[1]));
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
            }
    }
    return yMin;
}


vector<shared_ptr<Cell>> LimiterBJVertex::getStencilFor(const std::shared_ptr<Cell>& cell)
{
    vector<shared_ptr<Cell>> stencil = { cell }; // the stencil of limitation
    stencil.insert(stencil.end(), cell->neibCellsVertex.begin(), cell->neibCellsVertex.end());
    return stencil;
}


numvector<double, dimS> LimiterBJVertex::limitation(const std::vector<std::shared_ptr<Cell>>& stencil)
{
    int nCells = stencil.size();

    // load memory
    std::vector<numvector<double, dimPh>> uMean(nCells);    // mean values
    numvector<double, dimS> alphaNew;                       // result of limitation

    for (size_t k = 0; k < nCells; ++k)
        uMean[k] = solution.reconstruct(stencil[k]->number, stencil[k]->getCellCenter());

    // get minimum and maximum from cell averages
    
    numvector<double, dimPh> mI = 1e9;
    numvector<double, dimPh> MI = -1e9;

    for (int iSol = 0; iSol < dimPh; ++iSol)
    {
        auto q = std::minmax_element(uMean.begin(), uMean.end(), [&iSol](const numvector<double, dimPh>& a, const numvector<double, dimPh>& b){return a[iSol] < b[iSol]; });
        mI[iSol] = (*(q.first))[iSol];
        MI[iSol] = (*(q.second))[iSol];
    }

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
