#include "IndicatorEverywhere.h"

using namespace std;

vector<int> IndicatorEverywhere::checkDiscontinuities() const
{
    //vector<double> indicator(mesh.nCells);

    vector<int> troubledCells(mesh.nRealCells);
    int n = troubledCells.size();

//#pragma omp parallel for \
shared (n, troubledCells)

    for (size_t i = 0; i < n; ++i)
    {
        troubledCells[i] = i;

        for (int k = 0; k < dimPh; ++k)
            values[k*mesh.nRealCells + i] = 0.0;
    }

    return troubledCells;

}
