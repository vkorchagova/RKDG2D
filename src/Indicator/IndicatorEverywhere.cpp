#include "IndicatorEverywhere.h"

using namespace std;

vector<int> IndicatorEverywhere::checkDiscontinuities() const
{
    //vector<double> indicator(mesh.nCells);

    vector<int> troubledCells(mesh.nCells);

    for (size_t i = 0; i < troubledCells.size(); ++i)
        troubledCells[i] = i;

    return troubledCells;

}
