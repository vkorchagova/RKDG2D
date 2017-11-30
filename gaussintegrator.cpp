#include "gaussintegrator.h"

using namespace std;

GaussIntegrator::GaussIntegrator()
{
    double sqrt3 = 1.0/1.732050807568877;

    gPoints1D = { -sqrt3, sqrt3 };
    gPoints2D = { {-sqrt3, -sqrt3}, {sqrt3, -sqrt3}, {sqrt3, sqrt3}, {-sqrt3, sqrt3} };

    gWeights1D = { 1.0, 1.0 };
    gWeights2D = { 1.0, 1.0, 1.0, 1.0 };
}

// private

double GaussIntegrator::localToGlobal(double x, double a, double b)
{
    return 0.5*(b-a)*x + 0.5*(a+b);
}

numvector<double, 2> GaussIntegrator::localToGlobal(numvector<double, 2> coord, numvector<numvector<double,2>,4> nodes)
{
    double ax = nodes[0][0];
    double ay = nodes[0][1];
    double bx = nodes[1][0];
    double by = nodes[2][1];

    return { 0.5*(bx - ax)*coord[0] + 0.5*(bx + ax), 0.5*(by - ay)*coord[1] + 0.5*(by + ay) };
}

// ------------------------------------------ public

double GaussIntegrator::integrate(const function<double(double)>& f, double a, double b)
{
    int nGP = 2;

    double res = 0.0;
    double J = 0.5*(b-a);

    for (int i = 0; i < nGP; ++i)
    {
        res += gWeights1D[i]*f( localToGlobal(gPoints1D[i],a,b) );
    }

    return res*J;
} // for integrate 1D scalar

numvector<double,5> GaussIntegrator::integrate( const function<numvector<double,5>(double)>& f, double a, double b)
{
    int nGP = 2;

    numvector<double, 5> res (0.0);
    double J = 0.5*(b-a);

    for (int i = 0; i < nGP; ++i)
    {
        numvector<double, 5> resF = f( localToGlobal(gPoints1D[i],a,b) );

        for (int k = 0; k < 5; ++k)
        {
            res[k] += gWeights1D[i] * resF[k];
        }
    }

    return res * J;
} // for integrate 1D vector


double GaussIntegrator::integrate(const function<double(const numvector<double, 2>&)>& f, const numvector<numvector<double,2>,4>& nodes)
{
    int nGP = 4;

    double res = 0.0;
    double J = getArea(nodes)*0.25;

    for (int i = 0; i < nGP; ++i)
    {
        res += gWeights2D[i]*f(localToGlobal(gPoints2D[i],nodes));
    }

    return res*J;
} // for integrate 2D of scalar function

numvector<double,5> GaussIntegrator::integrate( const function<numvector<double,5>(const numvector<double, 2>&)>& f, const numvector<numvector<double,2>,4>& nodes)
{
    int nGP = 4;

    numvector<double,5> res = 0.0;
    double J = getArea(nodes)*0.25;

    for (int i = 0; i < nGP; ++i)
    {
        numvector<double,5> resF = f(localToGlobal(gPoints2D[i],nodes));
        for(int k=0;k<5;++k)
        {
            res[k] += gWeights2D[i] * resF[k];
        }
    }

    return res*J;
} // for integrate 2D of vector function

