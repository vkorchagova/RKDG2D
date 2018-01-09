#include "LimiterWENOS.h"

using namespace std;

void LimiterWENOS::limit(vector<numvector<double, 5 * nShapes>>& alpha)
{
    problem.setAlpha(alpha);
    vector<double> ind = indicator.checkDiscontinuities();

    // linear weights

    numvector<double, 3> gamma = { 0.001, 0.998, 0.001};

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

    for (size_t cell = 0; cell < alpha.size(); ++cell)
    {
        if (ind[cell] > 1.0)
        {
            // limit solution in X direction

            // find neighbours in x direction (we know the mesh is rectilinear)

            // boundary cells??????

            // norm of FF!!

            numvector<shared_ptr<Cell>, 3> cellsHor = { indicator.mesh.cells[cell-1], indicator.mesh.cells[cell], indicator.mesh.cells[cell+1] };

            // get mean values

            for (int i = 0; i < 3; ++i)
                uMean[i] = cellsHor[i]->h().x() / cellsHor[2]->h().x() * cellsHor[2]->reconstructSolution(cellsHor[i]->getCellCenter());

            // get coeffs for polynoms p

            p[0] = alpha[cell];

            for (int i = 0; i < 5; ++i)
                p[0][i*nShapes] += - uMean[0][i] + uMean[1][i];

            p[1] = alpha[cell];
            p[2] = alpha[cell];

            for (int i = 0; i < 5; ++i)
                p[2][i*nShapes] += - uMean[2][i] + uMean[1][i];

            // get smoothness indicators and nonlinear weights

            for (int i = 0; i < 3; ++i)
            {
                for (int j = 0; j < 5; ++j)
                {
                    beta[i][j] = sqr( cellsHor[i]->h().x() * alpha[cellsHor[i]->number][j*nShapes + 1]);
                    wTilde[i][j] = gamma[i] * 1.0 / sqr(beta[i][j] + 1e-6);
                }

                wSum = wTilde[0] + wTilde[1] + wTilde[2];

                for (int j = 0; j < 5; ++j)
                    w[i][j] = wTilde[i][j] / wSum[j];
            }

            // calc limited solution

            numvector<double, 5*nShapes> aLimited (0.0);

            for (int k = 0; k < 3; ++k)
                for (int j = 0; j < 5; ++j)
                    for (int q = 0; q < nShapes; ++q)
                        aLimited[j*nShapes + q] += w[k][j] * p[k][j*nShapes + q];

            alpha[cell] = aLimited;
        }
    }

    return;
}
