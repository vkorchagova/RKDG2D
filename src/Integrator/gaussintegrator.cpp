#include "gaussintegrator.h"

using namespace std;

GaussIntegrator::GaussIntegrator()
{
    double sqrt3 = 1.0/1.732050807568877;

    gPoints1D.set(-sqrt3, sqrt3);

    gPoints2D[0].set(-sqrt3, -sqrt3);
    gPoints2D[1].set( sqrt3, -sqrt3);
    gPoints2D[2].set(-sqrt3,  sqrt3);
    gPoints2D[3].set( sqrt3,  sqrt3);

    gWeights1D = { 1.0, 1.0 };
    gWeights2D = { 1.0, 1.0, 1.0, 1.0 };
}

// private

double GaussIntegrator::localToGlobal(double x, double a, double b)
{
    return 0.5*(b-a)*x + 0.5*(a+b);
}

Point& GaussIntegrator::localToGlobal(Point &point, Cell &cell)
{
    Point global (0.5 * cell.hx * point.x() + cell.center.x(), 0.5*cell.hy*point.y() + cell.center.y());

    return global;
}

// ------------------------------------------ public

/*
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
*/

numvector<double,5> GaussIntegrator::integrate( const std::function<numvector<double,5>(double)>& f, double a, double b)
{
    int nGP = 2;

    numvector<double, 5> res (0.0);
    double J = 0.5*(b-a);

    for (int i = 0; i < nGP; ++i)
    {
        numvector<double, 5> resF = f( localToGlobal(gPoints1D.x(),a,b) );

        for (int k = 0; k < 5; ++k)
        {
            res[k] += gWeights1D[i] * resF[k];
        }
    }

    return res * J;
} // for integrate 1D vector

/*
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
*/


numvector<double,5> GaussIntegrator::integrate( const std::function<numvector<double, 5>(const Point &)> &f, Cell& cell)
{
    int nGP = 4;

    numvector<double,5> res = 0.0;
    double J = cell.area*0.25;

    for (int i = 0; i < nGP; ++i)
    {
        numvector<double,5> resF = f(localToGlobal(gPoints2D[i],cell));
        for (int k = 0; k < 5; ++k)
        {
            res[k] += gWeights2D[i] * resF[k];
        }
    }

    return res*J;
} // for integrate 2D of vector function

