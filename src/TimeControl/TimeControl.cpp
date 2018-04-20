#include "TimeControl.h"
#include <algorithm>

using namespace std;

void TimeControl::getMassFlux()
{
    for (int i = 0; i < mesh.nCells; ++i)
        massFlux[i] += mesh.cells[i]->totalMassFlux();

}

void TimeControl::updateTimeStep()
{
    if (modifyTime)
    {
        double factCo = 0.0;
        double relTau = 0.0;
        vector<double> newTauLocal;
        newTauLocal.reserve(mesh.nCells);

        for (int i = 0; i < mesh.nCells; ++i)
        {
            factCo = 0.5 * tauOld * massFlux[i] / mesh.cells[i]->totalMass();
            relTau = CoNum / (factCo + 1e-6); //1e-6 is technical small number
            relTau = min(relTau, maxTauGrowth);
            tauNew = min(relTau * tauOld, maxTau);
            newTauLocal.push_back(tauNew);
        }

        tauNew = *min_element(newTauLocal.begin(),newTauLocal.end());
        tauOld = tauNew;
    }

    fill(massFlux.begin(), massFlux.end(), 0.0);
}
