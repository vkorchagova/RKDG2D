#include "LimiterFinDiff.h"

using namespace std;

void LimiterFinDiff::limit(vector<numvector<double, 5 * nShapes>>& alpha)
{
    problem.setAlpha(alpha);
    vector<int> troubledCells = indicator.checkDiscontinuities();

    for (int iCell : troubledCells)
    {
        for (int j = 0; j < 5; ++j)
        {
            alpha[iCell][j*nShapes + 1] = 0.0;
            alpha[iCell][j*nShapes + 2] = 0.0;
        }
    }
}
