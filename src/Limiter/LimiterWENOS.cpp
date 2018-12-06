#include "LimiterWENOS.h"
#include <omp.h>

using namespace std;

void LimiterWENOS::limit(vector<numvector<double, dimS>>& SOL)
{
    double ts = omp_get_wtime();
    int n = SOL.size();
    
    vector<numvector<double, dimS>> SOL_Corr(n);

#pragma omp parallel for \
shared(n, SOL_Corr, SOL) \
default(none)
    for (int i = 0; i < n; ++i)
		SOL_Corr[i] = SOL[i];
    
    // troubled cells

    vector<int> troubledCells;

    // linear weights

    vector<double> gamma;
    double g = 0.001;

    int nIter = 1;

    // smoothness indicators

    vector<numvector<double, dimPh>> beta;

    // nonlinear weights

    vector<numvector<double, dimPh>> w;
    vector<numvector<double, dimPh>> wTilde;
    numvector<double, dimPh> wSum;

    // mean values

    vector<numvector<double, dimPh>> uMean;

    // p polynoms

    vector<numvector<double, dimS>> p;


    troubledCells = indicator.checkDiscontinuities();

        // limit solution in troubled cells
#pragma omp parallel for \
shared(g, SOL, SOL_Corr, troubledCells) \
private(uMean, beta, gamma, w, wTilde, wSum, p) \
default(none)
        for (size_t i = 0; i < troubledCells.size(); ++i)
        //for (int iCell : troubledCells)
        {
			int iCell = troubledCells[i];

            shared_ptr<Cell> cell = M.cells[iCell];

            // construct list of cells: cell + neighbours

            vector<shared_ptr<Cell>> cells = { cell };
            cells.insert(cells.end(), cell->neibCells.begin(), cell->neibCells.end());

            int nCells = M.nCells;

            // get mean values of linear functions

            uMean.resize(nCells);

            for (size_t k = 0; k < nCells; ++k)
                uMean[k] = sln.reconstructSolution(cells[0]->getCellCenter());

            // get coeffs for polynoms p = a0 + a1 * (x - xc) + a2 * (y - yc)

            p.resize(nCells);

            for (size_t k = 0; k < nCells; ++k)
                p[k] = SOL[cells[k]->number] ;

            for (size_t k = 0; k < nCells; ++k)
                for (int i = 0; i < dimPh; ++i)
                    for (int j = 0; j < nShapes; ++j)
                        p[k][i*nShapes + j] *= cells[k]->offsetPhi[j];


            for (size_t k = 0; k < nCells; ++k)
                for (int i = 0; i < 5; ++i)
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
                for (int j = 0; j < 5; ++j)
                {
                    beta[k][j] =  cells[0]->getArea() * cells[k]->getArea()  * (sqr(p[k][j*nShapes + 1]) + sqr(p[k][j*nShapes + 2]));
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


            // project limited solution onto cell basis

            function<numvector<double, 5>(const Point& r)> foo = [=](const Point& r) \
            {
                numvector<double, 5> sum (0.0);

                for (int i = 0; i < 5; ++i)
                    for (size_t k = 0; k < nCells; ++k)
                        sum[i] += w[k][i] * (p[k][i*nShapes] + \
                                             p[k][i*nShapes + 1] * (r.x() - cells[k]->getCellCenter().x()) + \
                                             p[k][i*nShapes + 2] * (r.y() - cells[k]->getCellCenter().y()));


                return sum;
            };

            alphaNew[iCell] = cells[0]->projection(foo);
        }

    alpha = alphaNew;

#pragma omp parallel for \
shared(alpha, troubledCells, indicator, cout) \
default(none)
    for (size_t i = 0; i < troubledCells.size(); ++i)
    {
    //for (const shared_ptr<Cell> cell : indicator.mesh.cells)
        const shared_ptr<Cell> cell = indicator.mesh.cells[troubledCells[i]];
    
        for (const shared_ptr<Point> node : cell->nodes)
        {
            numvector<double, 5> res = cell->reconstructSolution(node);

            if (res[0] < 0 || res[4] < 0 || problem.getPressure(res) < 0)
            {
                cout << "negative values after limitation in cell #" << cell->number << endl;
                cout << "rho | rhoU | e = " << cell->reconstructSolution(node) << endl;
                cout << "p = " << problem.getPressure(cell->reconstructSolution(node)) << endl;

                for (int j = 0; j < 5; ++j)
                {
                    alpha[cell->number][j*nShapes + 1] = 0.0;
                    alpha[cell->number][j*nShapes + 2] = 0.0;
                }
            }
        }
        
    }

    //problem.setAlpha(alphaNew);
    
    
    double te = omp_get_wtime();
    
    //cout << "WENOS limiter " << te - ts << endl; 
}
