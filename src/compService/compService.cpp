#include "compService.h"
#include "defs.h"

using namespace std;

numvector<double, PhysDim> inverseRotate(const numvector<double, PhysDim>& sol, const Point& n)
{
    return { sol[0], n.x() * sol[1] - n.y() * sol[2],  n.y() * sol[1] + n.x() * sol[2], sol[3], sol[4] };
}

numvector<double, 5> rotate(const numvector<double,5> &sol, const Point& n)
{
    return { sol[0], n.x() * sol[1] + n.y() * sol[2], - n.y() * sol[1] + n.x() * sol[2], sol[3], sol[4] };
}


Point inverseRotate(const Point& v, const Point& n)
{
    return Point({ n.x() * v.x() - n.y() * v.y(),  n.y() * v.x() + n.x() * v.y() });
}

Point rotate(const Point& v, const Point& n)
{
    return Point({ n.x() * v.x() + n.y() * v.y(), - n.y() * v.x() + n.x() * v.y() });
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
numvector<double, PhysDim> integrate(const Cell& cell, const function<numvector<double, PhysDim>(const Point &)>& f)
{
    numvector<double, PhysDim> res (0.0);

    for (int i = 0; i < cell.nGP; ++i)
    {
         numvector<double,PhysDim> resF = f(cell.gPoints2D[i]);

        for (int k = 0; k < PhysDim; ++k)
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
	return distance(globalNumbers.begin(),find(globalNumbers.begin(), globalNumbers.end(), curNum));
}

int getPatchByName(std::vector<Patch>& patches, std::string& pName)
{
    for (int i = 0; i < patches.size(); ++i)
        if (patches[i].name == pName)
            return i;

    return -1;
}