#include "LimiterFinDiff.h"

using namespace std;

void LimiterFinDiff::limit(vector<numvector<double, 5 * nShapes>>& alpha)
{
    problem.setAlpha(alpha);
    vector<double> ind = indicator.checkDiscontinuities();

    for (size_t i = 0; i < alpha.size(); ++i)
    {
        if (ind[i] > 1.0)
        {
            for (int j = 0; j < 5; ++j)
            {
                alpha[i][j*nShapes + 1] = 0.0;
                alpha[i][j*nShapes + 2] = 0.0;
            }
        }
    }

    return;
}
