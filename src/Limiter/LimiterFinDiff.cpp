#include "LimiterFinDiff.h"

using namespace std;

vector<shared_ptr<Cell>> LimiterFinDiff::getStencilFor(const std::shared_ptr<Cell>& cell)
{
    vector<shared_ptr<Cell>> stencil = { cell }; 
    return stencil;
}

numvector<double, dimS> LimiterFinDiff::limitation(const std::vector<std::shared_ptr<Cell>>& stencil)
{
    numvector<double, dimS> res = solution.SOL[stencil[0]->number];
    
    for (int j = 0; j < dimS; ++j)
    {
        res[j*nShapes + 1] = 0.0;
        res[j*nShapes + 2] = 0.0;
    }

    return res;
}
