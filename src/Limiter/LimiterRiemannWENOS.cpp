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
            beta[k][j] =  cells[0]->getArea() * cells[k]->getArea() * (sqr(p[k][j*nShapes + 1]) + sqr(p[k][j*nShapes + 2]));
            wTilde[k][j] = gamma[k] * (1.0 / sqr(beta[k][j] + 1e-6));
        }

    wSum = {0.0, 0.0, 0.0, 0.0, 0.0};

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

    return cells[0]->projection(foo);
}

void LimiterRiemannWENOS::limit(vector<numvector<double, 5 * nShapes>>& alpha)
{
    problem.setAlpha(alpha);

    vector<int> troubledCells = indicator.checkDiscontinuities();

    // p polynoms: first() is px, second() is py

    vector<pair<numvector<double, 5 * nShapes>, numvector<double, 5 * nShapes>>> p;
    vector<numvector<double, 5 * nShapes>> px;
    vector<numvector<double, 5 * nShapes>> py;

    // limit solution in troubled cells

    for (int iCell : troubledCells)
    {
        shared_ptr<Cell> cell = indicator.mesh.cells[iCell];

        // construct list of cells: cell + neighbours

        vector<shared_ptr<Cell>> cells = { cell };
        cells.insert(cells.end(), cell->neibCells.begin(), cell->neibCells.end());

        int nCells = cells.size();

        // get px and py polynoms as Riemann invariants: first() is px, second() is py

        p.resize(nCells);
        px.resize(nCells);
        py.resize(nCells);

        for (int k = 0; k < nCells; ++k)
        {
            p[k] = cells[k]->getRiemannInvariants();
            px[k] = p[k].first;
            py[k] = p[k].second;
        }

        p[0].first = limitP(cells, px);
        p[0].second = limitP(cells, py);

        alpha[iCell] = cells[0]->reconstructCoefficients(p[0]);

        problem.setAlpha(alpha);

    }

}

