#include "RungeKutta.h"

using namespace std;

RungeKutta::RungeKutta(int o, Solver& s, Limiter& l, Time& t) : DDT(o, s, l, t)
{
    nSteps = order; // OK for RK order <= 4

    k.resize(nSteps);

    alpha.resize(nSteps);
    beta.resize(nSteps);

    for (int i = 0; i < nSteps; ++i)
        beta[i].resize(nSteps);

    setButcherTable();
}

void RungeKutta::setButcherTable()
{
    switch (order)
    {
        case 1:

            beta[0][0] = 1.0;
            break;

        case 2:

            alpha[0] = 0.0;
            alpha[1] = 1.0;

            //beta[0][0] = 0.5;
            //beta[1][0] = 0.0;
            //beta[1][1] = 1.0;

            beta[0][0] = 1.0;
            beta[1][0] = 0.5;
            beta[1][1] = 0.5;

            break;

        case 3:

            alpha[1] = 0.5;
            alpha[2] = 1.0;

            beta[0][0] = 0.5;
            beta[1][1] = 1.0;

            beta[2][0] = 0.1666666666666667;
            beta[2][1] = 0.6666666666666667;
            beta[2][2] = 0.1666666666666667;

            break;

        case 4:

            alpha[1] = 0.5;
            alpha[2] = 0.5;
            alpha[3] = 1.0;

            beta[0][0] = 0.5;
            //beta[1][0] = 0.0;
            beta[1][1] = 0.5;
            //beta[2][0] = 0.0;
            //beta[2][1] = 0.0;
            beta[2][2] = 1.0;

            beta[3][0] = 0.1666666666666667;
            beta[3][1] = 0.3333333333333333;
            beta[3][2] = 0.3333333333333333;
            beta[3][3] = 0.1666666666666667;

            break;

        default:
            cout << "Runge --- Kutta method of order " \
                 << order \
                 << "is not implemented. " \
                 << "Please change order to 1, 2, 3 or 4." \
                 << endl;
            exit(1);

    }
}

vector<numvector<double, 5 * nShapes>> RungeKutta::update\
         (vector<numvector<double, 5 * nShapes>> yOld, double tau)
{
    double tOld = time.runTime();

    vector<numvector<double, 5 * nShapes>> yNew = yOld;
    vector<numvector<double, 5 * nShapes>> lhs  = yOld;
    vector<numvector<double, 5 * nShapes>> lhsOld = solver.correctPrevIter(yOld);

    for (int i = 0; i < nSteps; ++i)
    {
        time.updateTime(tOld + alpha[i]*tau);
        k[i] = solver.assembleRHS(yNew);

        lhs = lhsOld;

        for (int j = 0; j <= i; ++j)
           lhs = lhs + k[j] * beta[i][j] * tau;

        yNew = solver.correctNonOrtho(lhs);

        limiter.limit(yNew);
    }
    
    //time.updateTime(tOld);

    return yNew;
}