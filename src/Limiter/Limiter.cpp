#include "Limiter.h"

using namespace std;

void Limiter::lastHope(std::vector<numvector<double, dimS> >& alpha)
{
    for (const shared_ptr<Cell>& cell : cells)
    {       
        for (const shared_ptr<Point>& node : cell->nodes)
        {
            numvector<double, dimPh> res = solution.reconstruct(cell->number, *node);

            if (res[0] < 0 || res[4] < 0 || physics.getPressure(res) < 0)
            {
//                cout << "negative values after limitation in cell #" << cell->number << endl;
//                cout << "rho | rhoU | e = " << cell->reconstructSolution(node) << endl;
//                cout << "p = " << problem.getPressure(cell->reconstructSolution(node)) << endl;

                for (int j = 0; j < dimPh; ++j)
                {
                    alpha[cell->number][j*nShapes + 1] = 0.0;
                    alpha[cell->number][j*nShapes + 2] = 0.0;
                }
            }
        }
    }
}