#include "compService.h"

using namespace std;

double integrate(const Cell& cell, const function<double(const Point &)>& f)
{
    double res = 0.0;

    for (int i = 0; i < cell.nGP; ++i)
    {
        double resF = f(cell.gPoints2D[i]);

        res += cell.gWeights2D[i] * resF * cell.J[i];

    }

    return res;

} // end integrate 2D of scalar function


numvector<double, 5> integrate(const Cell& cell, const function<numvector<double, 5>(const Point &)>& f)
{
    numvector<double, 5> res (0.0);

    for (int i = 0; i < cell.nGP; ++i)
    {
         numvector<double,5> resF = f(cell.gPoints2D[i]);

        for (int k = 0; k < 5; ++k)
        {
           res[k] += cell.gWeights2D[i] * resF[k] * cell.J[i];
        }
    }

    return res;

} // end integrate 2D of vector function
