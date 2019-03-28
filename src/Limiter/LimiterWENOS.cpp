#include "LimiterWENOS.h"
#include <omp.h>
#include <iostream>

using namespace std;

LimiterWENOS::LimiterWENOS
(
    const std::vector<std::shared_ptr<Cell>>& cells, 
    const Solution& sln,
    const Physics& phs
) : Limiter(cells, sln, phs) 
{
    

    //beta.reserve(maxPossibleStencilSize);
    //w.reserve(maxPossibleStencilSize);
    //wTilde.reserve(maxPossibleStencilSize);
    //uMean.reserve(maxPossibleStencilSize);
    //p.reserve(maxPossibleStencilSize);
}

void LimiterWENOS::limit(vector<numvector<double, dimS>>& alpha)
{
    int n = alpha.size();

    int maxPossibleStencilSize = 5;

    double t0 = MPI_Wtime();

    //troubledCells = indicator.checkDiscontinuities();
    troubledCells.resize(n);
#pragma omp parallel for \
 shared( n) \
 default(none)
    for (int i = 0; i < n; ++i)
        troubledCells[i] = i;  

    // stencil
    vector<shared_ptr<Cell>> stenc;
    vector<int> stencNumber;
    stenc.reserve(maxPossibleStencilSize);  
    stencNumber.reserve(maxPossibleStencilSize);

    /// Linear weights
    std::vector<double> gamma;

    /// smoothness indicators
    std::vector<numvector<double, dimPh>> beta;

    /// nonlinear weights
    std::vector<numvector<double, dimPh>> w;
    std::vector<numvector<double, dimPh>> wTilde;
    numvector<double, dimPh> wSum;

    /// mean values
    std::vector<numvector<double, dimPh>> uMean;

    /// p polynoms
    std::vector<numvector<double, dimS>> p;

    //beta.reserve(maxPossibleStencilSize);
    //w.reserve(maxPossibleStencilSize);
    //wTilde.reserve(maxPossibleStencilSize);
    //uMean.reserve(maxPossibleStencilSize);
    //p.reserve(maxPossibleStencilSize);


    // limit solution in troubled cells
#pragma omp parallel for shared(alpha, cout, cin) firstprivate(uMean, beta, gamma, w, wTilde, wSum, p, stenc, stencNumber) default(none)
    for (int iTroubled = 0; iTroubled < troubledCells.size(); ++iTroubled)
    //for (int iCell : troubledCells)
    {
        int iCell = troubledCells[iTroubled];

        shared_ptr<Cell> cell = cells[iCell];

        // construct list of cells: cell + neighbours

        stenc = { cell };
        stenc.insert(stenc.end(), cell->neibCells.begin(), cell->neibCells.end());

        int nCells = stenc.size();
        stencNumber.resize(nCells);
        for (size_t k = 0; k < nCells; ++k)
            stencNumber[k] = stenc[k]->number;

        // get mean values of linear functions

        uMean.resize(nCells);

        for (size_t k = 0; k < nCells; ++k)
            //uMean[k] = sln.reconstructSolution(cells[0]->getCellCenter());
            uMean[k] = solution.reconstruct(stencNumber[k], stenc[0]->getCellCenter());

        // get coeffs for polynoms p = a0 + a1 * (x - xc) + a2 * (y - yc)

        p.resize(nCells);

        for (size_t k = 0; k < nCells; ++k)
            p[k] = alpha[stencNumber[k]];
/*
cout << "after alpha" << endl;
cout << "cell #" << cell->number << ":\n";
        
for (size_t k = 0; k < nCells; ++k) \
                cout << "p #" << k << ": " << p[k] << endl;
*/
        for (size_t k = 0; k < nCells; ++k)
            for (int i = 0; i < dimPh; ++i)
                for (int j = 0; j < nShapes; ++j)
                    p[k][i*nShapes + j] *= solution.B.phiCoeffs[k][j];
/*
cout << "after phiCoeffs" << endl;
cout << "cell #" << cell->number << ":\n";
        
for (size_t k = 0; k < nCells; ++k) \
                cout << "p #" << k << ": " << p[k] << endl;
*/

        for (size_t k = 0; k < nCells; ++k)
            for (int i = 0; i < dimPh; ++i)
                p[k][i*nShapes] += - uMean[k][i] + uMean[0][i];
/*
if (iCell == 50)
{
cout << "after uMean" << endl;
            cout << "cell #" << cell->number << ":\n";
        
for (size_t k = 0; k < nCells; ++k) \
                cout << "p #" << k << ": " << p[k] << endl;
}*/

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
                beta[k][j] = 
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


        //cout << "beta:\n";
        //for (size_t k = 0; k < nCells; ++k)
        //{
        //    cout << "cell no = " << k << endl;
        //    cout << beta[k] << endl;
        //}

//        cout << "wtilde:\n";
//        for (size_t k = 0; k < nCells; ++k)
//        {
//            cout << "cell no = " << k << endl;
//            cout << wTilde[k] << endl;
//        }
/*if (iCell == 50)
{
        cout << "\tw:\n";
        for (size_t k = 0; k < nCells; ++k)
        {
            cout << "\tcell stenc no = " << k << "; ";
            cout << "\tcell no = " << stenc[k]->number << "; ";
            cout << "\tw[k] = " << w[k] << endl;
        }
    }*/



        function<numvector<double, dimPh> (const Point& r)> foo = [&](const Point& r) \
        {
            numvector<double, dimPh> sum (0.0);

            //for (size_t k = 0; k < nCells; ++k) \
                cout << "p #" << k << ": " << w[k] << endl;

            for (int i = 0; i < dimPh; ++i)
                for (size_t k = 0; k < nCells; ++k)
                    sum[i] += w[k][i] * (p[k][i*nShapes] + \
                                         p[k][i*nShapes + 1] * (r.x() - stenc[k]->getCellCenter().x()) + \
                                         p[k][i*nShapes + 2] * (r.y() - stenc[k]->getCellCenter().y()));


            return sum;
        };

        alphaNew[iCell] = solution.B.projection(foo, stencNumber[0]); //cells[0]->projection(foo);

        //if (iCell == 50)
         //   cout << alphaNew[iCell] << endl;
    }

    double t1 = MPI_Wtime();
    if (debug) logger << "\t\tLimiterWENOS.forTroubledCells: " << t1 - t0 << endl;

    alpha = alphaNew;

    t0 = MPI_Wtime();
    //lastHope(alpha);
    t1 = MPI_Wtime();
    if (debug) logger << "\t\tLimiterWENOS.lastHope: " << t1 - t0 << endl;
}
