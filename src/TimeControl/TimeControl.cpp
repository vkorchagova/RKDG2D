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

        for (const shared_ptr<Edge> edge : mesh.edges)
            edge->getMaxUL();

        for (const shared_ptr<Cell> cell : mesh.cells)
        {
            double uSum = 0.0;

            for (const shared_ptr<Edge> edge : cell->edges)
                uSum += edge->maxUL();

            factCo = 0.5 * tauOld * uSum / cell->getArea();
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
