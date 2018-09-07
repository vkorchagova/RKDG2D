#include "Adams.h"

using namespace std;

Adams::Adams(int o, Solver& s, Limiter& l, Time& t) : RungeKutta(o, s, l, t)
{
    nSteps = order; // OK for RK order <= 4
    b.resize(nSteps);

    setAdamsCoeffs();
}

void Adams::setAdamsCoeffs()
{
    switch (order)
    {
        case 1:

            b[0] = 1.0;
            break;

        case 2:

            b[0] =  1.5;
            b[1] = -0.5;

            break;

        case 3:

            b[0] =  23.0 / 12.0;
            b[1] = -16.0 / 12.0;
            b[2] =   5.0 / 12.0;

            break;

        case 4:

            b[0] =  55.0 / 24.0;
            b[1] = -59.0 / 24.0;
            b[2] =  37.0 / 24.0;
            b[3] = - 9.0 / 24.0;

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

vector<numvector<double, 5 * nShapes>> Adams::update\
         (vector<numvector<double, 5 * nShapes>> yOld, double tau)
{
    double tOld = time.runTime();
    vector<numvector<double, 5 * nShapes>> yNew = yOld;

    if (rhsHistory.size() < nSteps)
    {
        yNew = RungeKutta::update(yOld, tau);
        rhsHistory.push_front(solver.assembleRHS(yNew));
        time.updateTime(tOld + tau);
        return yNew;
    }

    vector<numvector<double, 5 * nShapes>> lhsOld = solver.correctPrevIter(yOld);
    vector<numvector<double, 5 * nShapes>> lhs  = lhsOld;

    

    for (int i = 0; i < nSteps; ++i)
        lhs = lhs + rhsHistory[i] * b[i] * tau;

    yNew = solver.correctNonOrtho(lhs);

    limiter.limit(yNew);

    rhsHistory.push_front(solver.assembleRHS(yNew));
    rhsHistory.pop_back();
    //time.updateTime(tOld + tau);

    

    return yNew;
}