#include "LimiterWENOS.h"

using namespace std;

void LimiterWENOS::limit(vector<numvector<double, 5 * nShapes>>& alpha)
{
    problem.setAlpha(alpha);
    vector<int> troubledCells = indicator.checkDiscontinuities();

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

    // p polynoms

    vector<numvector<double, 5 * nShapes>> p;


    // limit solution in troubled cells

    for (int iCell : troubledCells)
    {
        // find neighbours

        shared_ptr<Cell> cell = indicator.mesh.cells[iCell];

        vector<shared_ptr<Cell>> neibCellsX = cell->findNeighbourCellsX();
        vector<shared_ptr<Cell>> neibCellsY = cell->findNeighbourCellsY();

        vector<shared_ptr<Cell>> cells = { cell };


        cells.insert(cells.end(), neibCellsX.begin(), neibCellsX.end());

        cells.insert(cells.end(), neibCellsY.begin(), neibCellsY.end());

        int nCells = cells.size();

        // get mean values

        uMean.resize(nCells);

        for (size_t k = 0; k < nCells; ++k)
            uMean[k] = cells[k]->reconstructSolution(cells[0]->getCellCenter());

        // get coeffs for polynoms p

        for (size_t k = 0; k < nCells; ++k)
            p.emplace_back( alpha[cells[k]->number] );

        for (size_t k = 0; k < nCells; ++k)
            for (int i = 0; i < 5; ++i)
                for (int j = 0; j < nShapes; ++j)
                    p[k][i*nShapes + j] *= cells[k]->offsetPhi[j];

        for (size_t k = 1; k < nCells; ++k)
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
                beta[k][j] = cells[0]->h().x() * cells[0]->h().y() * (sqr(p[k][j*nShapes + 1]) + sqr(p[k][j*nShapes + 2]));
                wTilde[k][j] = gamma[k] * (1.0 / sqr(beta[k][j] + 1e-6));
            }

        wSum = {0.0, 0.0, 0.0, 0.0, 0.0};

        for (int j = 0; j < 5; ++j)
            for (size_t k = 0; k < nCells; ++k)
                wSum[j] += wTilde[k][j];

        for (size_t k = 0; k < nCells; ++k)
            for (int j = 0; j < 5; ++j)
                w[k][j] = wTilde[k][j] / wSum[j];

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

        alpha[iCell] = cells[0]->projection(foo);

        problem.setAlpha(alpha);

    }
}
