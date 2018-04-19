#include "TimeControl.h"
#include <algorithm>

using namespace std;

void TimeControl::updateTimeStep()
{
    if (modifyTime)
    {
        double factCo = 0.0;
        double relTau = 0.0;
        vector<double> newTauLocal;
        newTauLocal.reserve(mesh.nCells);

        for (const shared_ptr<Cell> cell : mesh.cells)
        {
            factCo = tauOld * cell->totalMassFlux() / cell->getArea();
            relTau = CoNum / (factCo + 1e-6); //1e-6 is technical small number
            relTau = min(relTau, maxTauGrowth);
            tauNew = min(relTau * tauOld, maxTau);
            newTauLocal.push_back(tauNew);
        }

        tauNew = distance(newTauLocal.begin(),min_element(newTauLocal.begin(),newTauLocal.end()));
        tauOld = tauNew;
    }
}
