#include "compService.h"
#include "defs.h"

using namespace std;

numvector<double, dimPh> inverseRotate(const numvector<double, dimPh>& sol, const Point& n)
{
    return { sol[0], n.x() * sol[1] - n.y() * sol[2],  n.y() * sol[1] + n.x() * sol[2], sol[3], sol[4] };
}

numvector<double, dimPh> rotate(const numvector<double, dimPh> &sol, const Point& n)
{
    return { sol[0], n.x() * sol[1] + n.y() * sol[2], - n.y() * sol[1] + n.x() * sol[2], sol[3], sol[4] };
}

numvector<double, dimGrad> rotate(const numvector<double, dimGrad> &gradSol, const Point& n)
{
    return { n.x() * gradSol[0] + n.y() * gradSol[1], -n.y() * gradSol[0] + n.x() * gradSol[1], \
                n.x() * ( gradSol[2] * n.x() + gradSol[4] * n.y()) + n.y() * ( gradSol[3] * n.x() + gradSol[5] * n.y() ), \
                -n.y() * ( gradSol[2] * n.x() + gradSol[4] * n.y()) + n.x() * ( gradSol[3] * n.x() + gradSol[5] * n.y()),\
                n.x() * ( gradSol[4] * n.x() - gradSol[2] * n.y()) + n.y() * ( gradSol[5] * n.x() - gradSol[3] * n.y() ), \
                -n.y() * ( gradSol[4] * n.x() - gradSol[2] * n.y()) + n.x() * ( gradSol[5] * n.x() - gradSol[3] * n.y()),\
                n.x() * gradSol[6] + n.y() * gradSol[7], -n.y() * gradSol[6] + n.x() * gradSol[7] };
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

double integrate(const Edge& edge, const function<double(const Point &)>& f)
{
    double res = 0.0;

    for (int i = 0; i < edge.nGP; ++i)
    {
        double resF = f(edge.gPoints[i]);

        res += edge.gWeights[i] * resF * edge.J;

    }

    return res;

} // end integrate 1D of scalar function


double integrate(const Edge& edge, const vector<double>& f)
{
    double res = 0.0;

    for (int i = 0; i < edge.nGP; ++i)
        res += edge.gWeights[i] * f[i] * edge.J;

    return res;

} // end integrate 1D of scalar function


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

numvector<double, dimPh> integrate(const Edge& edge, const function<numvector<double, dimPh>(const Point &)>& f)
{
    numvector<double, dimPh> res (0.0);

    for (int i = 0; i < edge.nGP; ++i)
    {
         numvector<double,dimPh> resF = f(edge.gPoints[i]);

        for (int k = 0; k < dimPh; ++k)
        {
           res[k] += edge.gWeights[i] * resF[k] * edge.J;
        }
    }

    return res;

} // end integrate 1D of vector function

numvector<double, dimPh> integrate(const Edge& edge, const vector<numvector<double, dimPh>>& f)
{
    numvector<double, dimPh> res (0.0);

    for (int i = 0; i < edge.nGP; ++i)
    {
        for (int k = 0; k < dimPh; ++k)
        {
           res[k] += edge.gWeights[i] * f[i][k] * edge.J;
        }
    }

    return res;

} // end integrate 1D of vector function



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

int getPatchByName(std::vector<ProcPatch>& patches, std::string& pName)
{
    for (int i = 0; i < patches.size(); ++i)
        if (patches[i].name == pName)
            return i;

    return -1;
}

///
/// Viscosity
///

numvector<double, dimGrad> outerProductArtificial(const numvector<double, dimPh>& sol, const Point& n)
{
	return
	{
		sol[0] * n[0],
		sol[0] * n[1],
		sol[1] * n[0],
		sol[1] * n[1],
		sol[2] * n[0],
		sol[2] * n[1],
		sol[4] * n[0],
		sol[4] * n[1]
	};
}
