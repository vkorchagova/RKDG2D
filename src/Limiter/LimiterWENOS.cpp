#include "LimiterWENOS.h"
#include <omp.h>
#include <iostream>

using namespace std;

LimiterWENOS::LimiterWENOS
(
    const Mesh& msh,
    Solution& sln,
    const Physics& phs,
    const Indicator& ind
) : Limiter(msh, sln, phs, ind) 
{
    // beta.reserve(maxPossibleStencilSize);
    // w.reserve(maxPossibleStencilSize);
    // wTilde.reserve(maxPossibleStencilSize);
    // uMean.reserve(maxPossibleStencilSize);
    // p.reserve(maxPossibleStencilSize);
    // stencNumber.reserve(maxPossibleStencilSize);
}

numvector<double, dimS> LimiterWENOS::limitP
(
    const vector<shared_ptr<Cell>>& stencil, 
    const vector<numvector<double, dimS>>& pInv
)
{
    int nCells = stencil.size(); 

    // load memory
    
    std::vector<double> gamma(nCells); // Linear weights
    std::vector<numvector<double, dimPh>> beta(nCells); // smoothness indicators
    std::vector<numvector<double, dimPh>> w(nCells); // nonlinear weights
    std::vector<numvector<double, dimPh>> wTilde(nCells);
    numvector<double, dimPh> wSum;
    std::vector<numvector<double, dimPh>> uMean(nCells); /// mean values
    

    /// p polynomas
    std::vector<numvector<double, dimS>> p = pInv;


    // get mean values of linear functions

    uMean.resize(nCells);

    for (size_t k = 0; k < nCells; ++k)
        //uMean[k] = sln.reconstructSolution(cells[0]->getCellCenter());
        uMean[k] = solution.reconstruct(stencil[k]->number, stencil[0]->getCellCenter(), p[k]);

    // get coeffs for polynoms p = a0 + a1 * (x - xc) + a2 * (y - yc)

    for (size_t k = 0; k < nCells; ++k)
        for (int i = 0; i < dimPh; ++i)
            for (int j = 0; j < nShapes; ++j)
                p[k][i*nShapes + j] *= solution.B.phiCoeffs[stencil[k]->number][j];
/*
cout << "after phiCoeffs" << endl;
cout << "cell #" << cell->number << ":\n";
    
for (size_t k = 0; k < nCells; ++k) \
            cout << "p #" << k << ": " << p[k] << endl;
*/

    for (size_t k = 0; k < nCells; ++k)
        for (int i = 0; i < dimPh; ++i)
            p[k][i*nShapes] += - uMean[k][i] + uMean[0][i];

    // get linear weights

    gamma.resize(nCells);

    gamma[0] = 1.0 - (nCells - 1) * g;

    for (size_t k = 1; k < nCells; ++k)
        gamma[k] = g;

    // get smoothness indicators and nonlinear weights

    beta.resize(nCells);
    wTilde.resize(nCells);
    w.resize(nCells);

    for (size_t k = 0; k < nCells; ++k)
        for (int j = 0; j < dimPh; ++j)
        {
            beta[k][j] = stencil[0]->getArea() * stencil[0]->getArea() *
                (sqr(p[k][j*nShapes + 1]) + sqr(p[k][j*nShapes + 2]));
            wTilde[k][j] = gamma[k] * (1.0 / sqr(beta[k][j] + 1e-6));
        }



    //wSum = {0.0, 0.0, 0.0, 0.0, 0.0};
    for (size_t k = 0; k < dimPh; ++k)
        wSum[k] = 0.0;


    for (int j = 0; j < dimPh; ++j)
        for (size_t k = 0; k < nCells; ++k)
            wSum[j] += wTilde[k][j];

//        //cout << wSum << endl;

    for (size_t k = 0; k < nCells; ++k)
        for (int j = 0; j < dimPh; ++j)
            w[k][j] = wTilde[k][j] / wSum[j];

//        cout << "----\n num tr cell = " << iCell << endl;

//        cout << "beta:\n";
//        for (size_t k = 0; k < nCells; ++k)
//        {
//            cout << "cell no = " << k << endl;
//            cout << beta[k] << endl;
//        }

//        cout << "wtilde:\n";
//        for (size_t k = 0; k < nCells; ++k)
//        {
//            cout << "cell no = " << k << endl;
//            cout << wTilde[k] << endl;
//        }

//        cout << "w:\n";
//        for (size_t k = 0; k < nCells; ++k)
//        {
//            cout << "cell no = " << k << endl;
//            cout << w[k] << endl;
//        }

    std::vector<numvector<double, dimS>> res(nCells);

    for (size_t k = 0; k < nCells; ++k)
        res[k] = numvector<double, dimS>(0.0);


    // for (int i = 0; i < dimPh; ++i)
    // {
    //     res[i*nShapes] =  p[0][i*nShapes];

    //     for (size_t k = 0; k < nCells; ++k)
    //     {
    //         res[i*nShapes + 1] += w[k][i] * p[k][i*nShapes + 1];
    //         res[i*nShapes + 2] += w[k][i] * p[k][i*nShapes + 2];
    //     }
    // }

    function<numvector<double, dimPh>(const Point& r)> foo = [=](const Point& r) \
    {
        numvector<double, dimPh> sum (0.0);

        for (int i = 0; i < dimPh; ++i)
            for (size_t k = 0; k < nCells; ++k)
                sum[i] += w[k][i] * (p[k][i*nShapes] + \
                                     p[k][i*nShapes + 1] * (r.x() - stencil[k]->getCellCenter().x()) + \
                                     p[k][i*nShapes + 2] * (r.y() - stencil[k]->getCellCenter().y()));


        return sum;
    };

    return solution.B.projection(foo, stencil[0]->number);
}


vector<shared_ptr<Cell>> LimiterWENOS::getStencilFor(const std::shared_ptr<Cell>& cell)
{
    vector<shared_ptr<Cell>> stencil = { cell }; 
    stencil.insert(stencil.end(), cell->neibCells.begin(), cell->neibCells.end());
    return stencil;
}


numvector<double, dimS> LimiterWENOS::limitation(const std::vector<std::shared_ptr<Cell>>& stencil)
{
    // see numbers of cells in the stencil

    int nCells = stencil.size();
    
    numvector<double, dimS> alphaNew; // result of limitation

    vector<numvector <double, dimS>> pInv(nCells);

    for (size_t k = 0; k < nCells; ++k)
        pInv[k] = solution.SOL[stencil[k]->number];

    //vector<numvector <double, dimS>> pLim = 

    

    alphaNew = limitP(stencil, pInv);; 

    return alphaNew;
}
