#include "LimiterWENOS.h"

using namespace std;

void LimiterWENOS::limit(vector<numvector<double, 5 * nShapes>>& alpha)
{
    problem.setAlpha(alpha);
    vector<int> troubledCells = indicator.checkDiscontinuities();

    // linear weights

    numvector<double, 3> gamma = { 0.001, 0.998, 0.001 };

    // smoothness indicators

    numvector<numvector<double, 5>, 3> beta;

    // nonlinear weights

    numvector<numvector<double, 5>, 3> w;
    numvector<numvector<double, 5>, 3> wTilde;
    numvector<double, 5> wSum;

    // mean values

    numvector<numvector<double, 5>, 3> uMean;

    // p polynoms

    numvector<numvector<double, 5 * nShapes>, 3> p;


    // limit solution in troubled cells

    for (int iCell : troubledCells)
    {
        // limit solution in X direction

        // find neighbours in x direction (we know the mesh is rectilinear)

        numvector<shared_ptr<Cell>, 3> cellsHor = { indicator.mesh.cells[iCell-1], \
                                                    indicator.mesh.cells[iCell], \
                                                    indicator.mesh.cells[iCell+1] };

        // get mean values

        for (int k = 0; k < 3; ++k)
            uMean[k] = cellsHor[k]->reconstructSolution(cellsHor[1]->getCellCenter());

        // get coeffs for polynoms p

        p[0] = alpha[iCell-1];
        p[1] = alpha[iCell];
        p[2] = alpha[iCell+1];

        for (int k = 0; k < 3; ++k)
            for (int i = 0; i < 5; ++i)
                for (int j = 0; j < nShapes; ++j)
                    p[k][i*nShapes + j] *= cellsHor[k]->offsetPhi[j];

        for (int i = 0; i < 5; ++i)
            p[0][i*nShapes] += - uMean[0][i] + uMean[1][i];

        for (int i = 0; i < 5; ++i)
            p[2][i*nShapes] += - uMean[2][i] + uMean[1][i];

        // get smoothness indicators and nonlinear weights

        for (int k = 0; k < 3; ++k)
        {
            for (int j = 0; j < 5; ++j)
            {
                beta[k][j] = sqr( cellsHor[k]->h().x() * p[k][j*nShapes + 1]);
                wTilde[k][j] = gamma[k] * (1.0 / sqr(beta[k][j] + 1e-6));
            }
        }

        wSum = wTilde[0] + wTilde[1] + wTilde[2];

        for (int k = 0; k < 3; ++k)
        {
            for (int j = 0; j < 5; ++j)
                w[k][j] = wTilde[k][j] / wSum[j];
        }

        // project limited solution onto cell basis

        function<numvector<double, 5>(const Point& r)> foo = [=](const Point& r) \
        {
            numvector<double, 5> sum (0.0);

            for (int i = 0; i < 5; ++i)
                for (int k = 0; k < 3; ++k)
                    sum[i] += w[k][i] * (p[k][i*nShapes] + \
                                         p[k][i*nShapes + 1] * (r.x() - cellsHor[k]->getCellCenter().x()) + \
                                         p[k][i*nShapes + 2] * (r.y() - cellsHor[k]->getCellCenter().y()));


            return sum;
        };

        alpha[iCell] = cellsHor[1]->projection(foo);

        problem.setAlpha(alpha);
    }
}
