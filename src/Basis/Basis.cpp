#include "Basis.h"
#include <vector>
#include <compService.h>

using namespace std;

Basis::Basis(const std::vector<shared_ptr<Cell>>& cls) : cells(cls)
{
    initBasisFunctions();
    initGramian();
}

void Basis::initBasisFunctions()
{
    for (const shared_ptr<Cell> cell : cells)
    {
        numvector<double,nShapes> curCoeffs =
        {
            1.0 / sqrt(cell->getArea()),
            2.0 * sqrt(3) / cell->getArea(),
            2.0 * sqrt(3) / cell->getArea()
        };

        if (nShapes == 6)
        {
            //curCoeffs.push_back(TODO)
        }

        phiCoeffs.push_back(curCoeffs);
    }

    phi.reserve(nShapes);
    gradPhi.reserve(nShapes);

    phi.emplace_back([&](int iCell, const Point& r){ return phiCoeffs[iCell][0]; });
    phi.emplace_back([&](int iCell, const Point& r){ return (r.x() - cells[iCell]->center.x()) * phiCoeffs[iCell][1]; });
    phi.emplace_back([&](int iCell, const Point& r){ return (r.y() - cells[iCell]->center.y()) * phiCoeffs[iCell][2]; });

    gradPhi.emplace_back([&](int iCell, const Point& r)->Point { return Point({ 0.0                , 0.0 }); });
    gradPhi.emplace_back([&](int iCell, const Point& r)->Point { return Point({ phiCoeffs[iCell][1], 0.0 }); });
    gradPhi.emplace_back([&](int iCell, const Point& r)->Point { return Point({ 0.0                , phiCoeffs[iCell][2] }); });

    if (nShapes == 6)
    {
        // CHECK FOR UNSTRUCTURED MESHES!!!
        //just for compilation, not used!!!
        double hx = 1.0;
        double hy = 1.0;

/*

        phi.emplace_back([&](const Point& r)
            { return sqrt(180.0 / pow(hx, 5) / hy) * ( sqr(r.x() - center.x()) - sqr(hx) / 12.0); });
        phi.emplace_back([&](const Point& r)
            { return sqrt(144.0 / pow(hx, 3) / pow(hy, 3)) * (r.x() - center.x()) * (r.y() - center.y()); });
        phi.emplace_back([&](const Point& r)
            { return sqrt(180.0 / pow(hy, 5) / hx) * ( sqr(r.y() - center.y()) - sqr(hy) / 12.0); });

        gradPhi.emplace_back([&](const Point& r)->Point
            { return Point({ sqrt(180.0 / pow(hx, 5) / hy) * 2 * (r.x() - center.x()) , 0.0}); });
        gradPhi.emplace_back([&](const Point& r)->Point
            { return Point({ sqrt(144.0 / pow(hx, 3) / pow(hy, 3)) * (r.y() - center.y()), \
                             sqrt(144.0 / pow(hx, 3) / pow(hy, 3)) * (r.x() - center.x()) }); });
        gradPhi.emplace_back([&](const Point& r)->Point
            { return Point({ 0.0, sqrt(180.0 / pow(hy, 5) / hx) * 2.0 * (r.y() - center.y()) }); });
*/
    }
}


void Basis::initGramian()
{
    gramian.reserve(cells.size());

    std::function<double(const Point&)> f;
    numvector<double, nShapes> curGram;


    for (int iCell = 0; iCell < cells.size(); ++iCell)
    {
        for (int i = 1; i < nShapes; ++i)
        {
            for (int j = i; j < nShapes; ++j)
            {

                f = [&](const Point& p) {  return phi[i](iCell,p) * phi[j](iCell,p); };
                curGram[(i - 1)*nShapes + j] = integrate(*cells[iCell], f);
            }
        }

        gramian.push_back(curGram);
    }
}

