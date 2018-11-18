#include "LimiterFinDiff.h"

using namespace std;

void LimiterFinDiff::limit(vector<numvector<double, 5 * nShapes>>& alpha)
{
    problem.setAlpha(alpha);
    vector<int> troubledCells = indicator.checkDiscontinuities();
    int n = troubledCells.size();

    //for (int iCell : troubledCells)
//#pragma omp parallel for \
    //for (vector<int>::iterator iter = troubledCells.begin(); iter != troubledCells.end(); ++iter)
    for (size_t i = 0; i < n; ++i)
    {        
        int iCell = troubledCells[i];//*iter;
        
        for (int j = 0; j < 5; ++j)
        {
            alpha[iCell][j*nShapes + 1] = 0.0;
            alpha[iCell][j*nShapes + 2] = 0.0;
        }
    }
}
