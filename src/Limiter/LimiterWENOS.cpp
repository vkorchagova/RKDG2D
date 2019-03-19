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
    int maxPossibleStencilSize = 4;

    beta.reserve(maxPossibleStencilSize);
    w.reserve(maxPossibleStencilSize);
    wTilde.reserve(maxPossibleStencilSize);
    uMean.reserve(maxPossibleStencilSize);
    p.reserve(maxPossibleStencilSize);
}

void LimiterWENOS::limit(vector<numvector<double, dimS>>& alpha)
{
    double ts = omp_get_wtime();
    int n = alpha.size();

    double t0 = MPI_Wtime();

    //troubledCells = indicator.checkDiscontinuities();
    troubledCells.resize(n);
#pragma omp parallel for \
 shared(troubledCells, n) \
 default(none)
    for (int i = 0; i < n; ++i)
        troubledCells[i] = i;  

    // stencil
    vector<shared_ptr<Cell>> stenc;
    stenc.reserve(5);  

        // limit solution in troubled cells
#pragma omp parallel for \
shared(alpha, alphaNew, troubledCells, cout) \
private(uMean, beta, gamma, w, wTilde, wSum, p, stenc) \
default(none)
    for (int iTroubled = 0; iTroubled < troubledCells.size(); ++iTroubled)
    //for (int iCell : troubledCells)
    {
        int iCell = troubledCells[iTroubled];

        shared_ptr<Cell> cell = cells[iCell];

        // construct list of cells: cell + neighbours

        stenc = { cell };
        stenc.insert(stenc.end(), cell->neibCells.begin(), cell->neibCells.end());

        int nCells = stenc.size();

        // get mean values of linear functions

        uMean.resize(nCells);

        for (size_t k = 0; k < nCells; ++k)
            //uMean[k] = sln.reconstructSolution(cells[0]->getCellCenter());
            uMean[k] = solution.reconstruct(stenc[k]->number, stenc[0]->getCellCenter());

        // get coeffs for polynoms p = a0 + a1 * (x - xc) + a2 * (y - yc)

        p.resize(nCells);

        for (size_t k = 0; k < nCells; ++k)
            p[k] = alpha[stenc[k]->number];

//cout << "after alpha" << endl;
//        for (size_t k = 0; k < nCells; ++k)
///                cout << "p #" << k << ": " << p[k] << endl;

        for (size_t k = 0; k < nCells; ++k)
            for (int i = 0; i < dimPh; ++i)
                for (int j = 0; j < nShapes; ++j)
                    p[k][i*nShapes + j] *= solution.B.phiCoeffs[k][j];

//cout << "after phiCoeffs" << endl;
//        for (size_t k = 0; k < nCells; ++k)
//                cout << "p #" << k << ": " << p[k] << endl;


        for (size_t k = 0; k < nCells; ++k)
            for (int i = 0; i < 5; ++i)
                p[k][i*nShapes] += - uMean[k][i] + uMean[0][i];

//cout << "after uMean" << endl;
//for (size_t k = 0; k < nCells; ++k)
//                cout << "p #" << k << ": " << p[k] << endl;



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
            for (int j = 0; j < 5; ++j)
            {
                beta[k][j] =  stenc[0]->getArea() * stenc[k]->getArea()  * (sqr(p[k][j*nShapes + 1]) + sqr(p[k][j*nShapes + 2]));
                wTilde[k][j] = gamma[k] * (1.0 / sqr(beta[k][j] + 1e-6));
            }



        //wSum = {0.0, 0.0, 0.0, 0.0, 0.0};
        wSum[0] = 0.0;
        wSum[1] = 0.0;
        wSum[2] = 0.0;
        wSum[3] = 0.0;
        wSum[4] = 0.0;


        for (int j = 0; j < 5; ++j)
            for (size_t k = 0; k < nCells; ++k)
                wSum[j] += wTilde[k][j];

//        //cout << wSum << endl;

        for (size_t k = 0; k < nCells; ++k)
            for (int j = 0; j < 5; ++j)
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

        //cout << "w:\n";
        //for (size_t k = 0; k < nCells; ++k)
        //{
        //    cout << "cell stenc no = " << k << ";";
        //    cout << "cell no = " << stenc[k]->number << ";";
        //    cout << "w[k] = " << w[k] << endl;
        //}

        // project limited solution onto cell B

        for (size_t k = 0; k < nCells; ++k)
                cout << "out p #" << k << ": " << p[k] << endl;

        function<numvector<double, 5>(const Point& r)> foo = [&](const Point& r) \
        {
            numvector<double, 5> sum (0.0);

            for (size_t k = 0; k < nCells; ++k)
                cout << "p #" << k << ": " << p[k] << endl;

            for (int i = 0; i < 5; ++i)
                for (size_t k = 0; k < nCells; ++k)
                    sum[i] += w[k][i] * (p[k][i*nShapes] + \
                                         p[k][i*nShapes + 1] * (r.x() - stenc[k]->getCellCenter().x()) + \
                                         p[k][i*nShapes + 2] * (r.y() - stenc[k]->getCellCenter().y()));


            return sum;
        };


        alphaNew[iCell] = solution.B.projection(foo,stenc[0]->number); //cells[0]->projection(foo);
    }
    double t1 = MPI_Wtime();
    if (debug) logger << "\t\tLimiterWENOS.forTroubledCells: " << t1 - t0 << endl;

    alpha = alphaNew;

    t0 = MPI_Wtime();
    lastHope(alpha);
    t1 = MPI_Wtime();
    if (debug) logger << "\t\tLimiterWENOS.lastHope: " << t1 - t0 << endl;
}
