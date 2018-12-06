#include "compService.h"

using namespace std;

numvector<double, dimPh> inverseRotate(const numvector<double, dimPh>& sol, const Point& n)
{
    return { sol[0], n.x() * sol[1] - n.y() * sol[2],  n.y() * sol[1] + n.x() * sol[2], sol[3], sol[4] };
}
numvector<double, 5> rotate(const numvector<double,5> &sol, const Point& n)
{
    return { sol[0], n.x() * sol[1] + n.y() * sol[2], - n.y() * sol[1] + n.x() * sol[2], sol[3], sol[4] };
}

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
numvector<double, dimPh> integrate(const Cell& cell, const function<numvector<double, dimPh>(const Point &)>& f)
{
    numvector<double, dimPh> res (0.0);

    for (int i = 0; i < cell.nGP; ++i)
    {
         numvector<double,dimPh> resF = f(cell.gPoints2D[i]);

        for (int k = 0; k < dimPh; ++k)
        {
           res[k] += cell.gWeights2D[i] * resF[k] * cell.J[i];
        }
    }

    return res;

} // end integrate 2D of vector function


///
/// MPI operations
///

int localNumber(std::vector<int>& globalNumbers, int curNum)
{
	return *find(globalNumbers.begin(), globalNumbers.end(), curNum);
}
