#include "Limiter.h"
#include <omp.h>

using namespace std;


Limiter::Limiter(
    const std::vector<std::shared_ptr<Cell>>& cells, 
    const Solution& sln,
    const Physics& phs) 
: cells(cells), solution(sln), physics(phs)
{
    int n = cells.size();

    alphaNew.resize(n);
    troubledCells.reserve(n);
}


void Limiter::lastHope(std::vector<numvector<double, dimS> >& alpha)
{
#pragma omp parallel for shared(alpha)
    for (int i = 0; i < cells.size(); ++i)
	//for (const shared_ptr<Cell>& cell : cells)
    {       
		const shared_ptr<Cell> cell = cells[i];
		
		for (const shared_ptr<Point>& node : cell->nodes)
        {
            numvector<double, dimPh> res = solution.reconstruct(cell->number, *node);
            double pres = physics.getPressure(res);

            if (res[0] < 0 || res[4] < 0 || pres < 0)
            {

                //cout << "\tnegative values after limitation in cell #" << cell->number << endl;
                //cout << "\t\trho = " << res[RHO] << endl;
                //cout << "\t\trhoU = " << res[RHOU] << endl;
                //cout << "\t\trhoV = " << res[RHOV] << endl;
                //cout << "\t\te = " << res[E] << endl;
                //cout << "\t\tp = " << pres << endl;

                for (int j = 0; j < dimPh; ++j)
                {
                    alpha[cell->number][j*nShapes + 1] = 0.0;
                    alpha[cell->number][j*nShapes + 2] = 0.0;
                }
            }
        }
    }
}