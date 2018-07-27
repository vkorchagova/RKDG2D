#include "LimiterRiemannWENOS.h"

using namespace std;

numvector<double, 5 * nShapes> limitP(const vector<shared_ptr<Cell>>& cells, const vector<numvector<double, 5 * nShapes>>& alpha)
{   // just for limiting of px and py

    int nCells = cells.size();
    vector<numvector<double, 5 * nShapes>> p = alpha;

    // linear weights

    vector<double> gamma;
    double g = 0.001;

    // smoothness indicators

    vector<numvector<double, 5>> beta;

    // nonlinear weights

    vector<numvector<double, 5>> w;
    vector<numvector<double, 5>> wTilde;
    numvector<double, 5> wSum;

    // mean values

    vector<numvector<double, 5>> uMean;

    // get mean values of Riemann invariants on the troubled cell

    uMean.resize(nCells);

    for (size_t k = 0; k < nCells; ++k)
        uMean[k] = cells[k]->reconstructSolution(cells[0]->getCellCenter(), p[k]);


    // get coeffs for polynoms p = a0 + a1 * (x - xc) + a2 * (y - yc)

    for (size_t k = 0; k < nCells; ++k)
        for (int i = 0; i < 5; ++i)
            for (int j = 0; j < nShapes; ++j)
                p[k][i * nShapes + j] *= cells[k]->offsetPhi[j];

    // move polynoms

    for (size_t k = 0; k < nCells; ++k)
        for (int i = 0; i < 5; ++i)
            p[k][i * nShapes] += - uMean[k][i] + uMean[0][i];

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
            beta[k][j] =  cells[0]->getArea() * (sqr(p[k][j*nShapes + 1]) + sqr(p[k][j*nShapes + 2]));
            wTilde[k][j] = gamma[k] * (1.0 / sqr(beta[k][j] + 1e-6));
        }

    wSum = {0.0, 0.0, 0.0, 0.0, 0.0};

    for (int j = 0; j < 5; ++j)
        for (size_t k = 0; k < nCells; ++k)
            wSum[j] += wTilde[k][j];

//        cout << wSum << endl;

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

//    function<numvector<double, 5>(const Point& r)> foo = [=](const Point& r) \
//    {
//        numvector<double, 5> sum (0.0);

//        for (int i = 0; i < 5; ++i)
//            for (size_t k = 0; k < nCells; ++k)
//                sum[i] += w[k][i] * (p[k][i*nShapes] + \
//                                     p[k][i*nShapes + 1] * (r.x() - cells[k]->getCellCenter().x()) + \
//                                     p[k][i*nShapes + 2] * (r.y() - cells[k]->getCellCenter().y()));


//        return sum;
//    };

//    return cells[0]->projection(foo);

    numvector<double, 5*nShapes> res (0.0);

    for (int i = 0; i < 5; ++i)
    {
        res[i*nShapes] =  p[0][i*nShapes];

        for (size_t k = 0; k < nCells; ++k)
        {
            //res[i*nShapes]     += w[k][i] * p[k][i*nShapes];
            res[i*nShapes + 1] += w[k][i] * p[k][i*nShapes + 1];
            res[i*nShapes + 2] += w[k][i] * p[k][i*nShapes + 2];
        }
    }

    return res;
}

void LimiterRiemannWENOS::limit(vector<numvector<double, 5 * nShapes>>& alpha)
{
    problem.setAlpha(alpha);

    vector<int> troubledCells = indicator.checkDiscontinuities();

    // p polynoms

    vector<numvector<double, 5 * nShapes>> p;
    vector<numvector<double, 5 * nShapes>> pInv;

    // limit solution in troubled cells

    for (int iCell : troubledCells)
    {
        shared_ptr<Cell> cell = indicator.mesh.cells[iCell];

        // construct list of cells: cell + neighbours

        vector<shared_ptr<Cell>> cells = { cell };
        cells.insert(cells.end(), cell->neibCells.begin(), cell->neibCells.end());

        int nCells = cells.size();

        // get p polynoms for each cell

        p.resize(nCells);
        pInv.resize(nCells);
        vector<numvector<double, 5 * nShapes>> pNew;

        for (int k = 0; k < nCells; ++k)
            p[k] = alpha[cells[k]->number];

        // get pNew polynoms using each edge normal
        for (const shared_ptr<Edge> e : cells[0]->edges)
        {
            if (e->neibCells.size() == 2)
            {
                // correct normal direction: outside related to troubled cell
                Point n = (( *e->nodes[0] - cell->getCellCenter()) * e->n > 0.0) ? e->n : Point(-e->n);

                // get Riemann invariants along this direction
                for (size_t k = 0; k < cells.size(); ++k)
                    pInv[k] = cells[k]->getRiemannInvariants(n);

                // limit Riemann invariants
                numvector <double, 5*nShapes> pLim = limitP(cells, pInv);

                // project limited solution to conservative variables
                pNew.push_back( cell->reconstructCoefficients(pLim, n) );
            }
        }

        double sumArea = 0.0;

        for (int i = 1; i < nCells; ++i)
            sumArea += cells[i]->getArea();

        // construct limited solution
        numvector<double, 5 * nShapes> res (0.0);
        for (int i = 0; i < nCells-1; ++i)
            res += pNew[i] * (cells[i+1]->getArea() / sumArea);


        function<numvector<double, 5>(const Point& r)> foo = [=](const Point& r) \
        {
            numvector<double, 5> sum (0.0);

            for (int i = 0; i < 5; ++i)
                sum[i] += (res[i*nShapes] + \
                           res[i*nShapes + 1] * (r.x() - cell->getCellCenter().x()) + \
                           res[i*nShapes + 2] * (r.y() - cell->getCellCenter().y()));


            return sum;
        };

        numvector<double, 5*nShapes> newAlpha = cell->projection(foo);
        alpha[cell->number] = newAlpha;
        problem.setAlpha(alpha);

    }

    for (const shared_ptr<Cell> cell : indicator.mesh.cells)
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

    problem.setAlpha(alpha);
}


