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

#pragma omp parallel for \
shared (mesh) \
default (none)
        for (size_t i = 0; i < mesh.edges.size(); ++i)
        {
            const shared_ptr<Edge> edge = mesh.edges[i];
            edge->getMaxUL();
        }

//        for (const shared_ptr<Edge> edge : mesh.edges)
//            edge->getMaxUL();

#pragma omp parallel for \
shared (mesh, newTauLocal) \
private (factCo, relTau) \
default (none)
        //for (const shared_ptr<Cell> cell : mesh.cells)
        for (size_t i = 0; i < mesh.cells.size(); ++i)
        {
            const shared_ptr<Cell> cell = mesh.cells[i];

            double uSum = 0.0;

            for (const shared_ptr<Edge> edge : cell->edges)
                uSum += edge->maxUL();

            factCo = tauOld * uSum / cell->getArea();
            relTau = CoNum / (factCo + 1e-6); //1e-6 is technical small number
            relTau = min(relTau, maxTauGrowth);
            tauNew = min(relTau * tauOld, maxTau);
            newTauLocal.push_back(tauNew);
        }

        tauNew = *min_element(newTauLocal.begin(),newTauLocal.end());

        tauOld = tauNew;
        fill(massFlux.begin(), massFlux.end(), 0.0);
    }


}
