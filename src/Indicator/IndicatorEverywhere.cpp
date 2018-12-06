#include "IndicatorEverywhere.h"

using namespace std;

vector<int> IndicatorEverywhere::checkDiscontinuities() const
{
    //vector<double> indicator(mesh.nCells);

    vector<int> troubledCells(mesh.nRealCells);
    int n = troubledCells.size();

    for (size_t i = 0; i < n; ++i)
        troubledCells[i] = i;

    return troubledCells;

}
