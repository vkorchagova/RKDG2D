#include <iostream>
#include <omp.h>

#include "Physics.h"
#include "Params.h"

using namespace std;

// ------------------ Constructors & Destructors ----------------

Physics::Physics()
{
    //cpcv = 1.4; // default value
} // end constructor 

Physics::~Physics()
{

}


// ------------------ Private class methods --------------------

// ------------------ Public class methods --------------------



double Physics::getPressure(const numvector<double, dimPh>& sol) const
{
    double rhoV2 = sqr(sol[1]) + sqr(sol[2]) + sqr(sol[3]);

    double p = (cpcv - 1.0)*(sol[4] - 0.5*rhoV2 / sol[0]);

    return p;
} // end getPressure

double Physics::c(const numvector<double, dimPh>& sol) const
{
    double c2 = cpcv * getPressure(sol) / sol[0];

    if (c2 < 0)
        cout << "sound speed < 0 =( !!!!" << endl;

    return sqrt( c2 );
} // end c for cell

numvector<double, dimPh> Physics::fluxF(const numvector<double, dimPh>& sol) const
{
    double u = sol[1] / sol[0];
    double p = getPressure(sol);

    return { sol[1], u*sol[1] + p, u*sol[2], u*sol[3], (sol[4] + p)*u };
} // end fluxF

numvector<double, dimPh> Physics::fluxG(const numvector<double, dimPh>& sol) const
{
    double v = sol[2] / sol[0];
    double p = getPressure(sol);

    return { sol[2], v*sol[1], v*sol[2] + p, v*sol[3], (sol[4] + p)*v };
} // end fluxG

numvector<numvector<double, dimPh>, dimPh> Physics::getLx(const numvector<double, dimPh>& sol) const
{
	double cS = c(sol);
	double rho = sol[0];
	double u = sol[1] / rho;
	double v = sol[2] / rho;

	double B1 = (cpcv - 1.0) / sqr(cS);
	double B2 = 0.5 * B1 * (sqr(u) + sqr(v));

	numvector<numvector<double, dimPh>, dimPh> res;
	res =
	{
		// Lx
		{ 0.5 * (B2 + u / cS), -0.5 * (B1 * u + 1.0 / cS), -0.5 * B1 * v,  0.0, 0.5 * B1 },
		{ -v,                   0.0,						1.0,		   0.0,      0.0 },
		{ 0.0,                  0.0,						0.0,		   1.0,      0.0 },
		{ 1.0 - B2,             B1 * u,						B1 * v,		   0.0,      -B1 },
		{ 0.5 * (B2 - u / cS), -0.5 * (B1 * u - 1.0 / cS), -0.5 * B1 * v,  0.0, 0.5 * B1 }
	};
	return res;
}
numvector<numvector<double, dimPh>, dimPh> Physics::getLy(const numvector<double, dimPh>& sol) const
{
	double cS = c(sol);
	double rho = sol[0];
	double u = sol[1] / rho;
	double v = sol[2] / rho;

	double B1 = (cpcv - 1.0) / sqr(cS);
	double B2 = 0.5 * B1 * (sqr(u) + sqr(v));

	numvector<numvector<double, dimPh>, dimPh> res;
	res =
	{
		// Ly
		{ 0.5 * (B2 + v / cS), -0.5 * (B1 * u), -0.5 * (B1 * v + 1.0 / cS),  0.0, 0.5 * B1 },
		{ u,				   -1.0,             0.0,						 0.0,      0.0 },
		{ 0.0,					0.0,             0.0,						 1.0,      0.0 },
		{ 1.0 - B2,				B1 * u,          B1 * v,					 0.0,      -B1 },
		{ 0.5 * (B2 - v / cS), -0.5 * (B1 * u), -0.5 * (B1 * v - 1.0 / cS),  0.0, 0.5 * B1 }
	};

	return res;
}
numvector<numvector<double, dimPh>, dimPh> Physics::getRx(const numvector<double, dimPh>& sol) const
{
    double cS = c(sol);
    double rho = sol[0];
    double u = sol[1] / rho;
    double v = sol[2] / rho;
    double p = getPressure(sol);

    double H = (sol[4] + p) / rho;
    double magU2 = 0.5 * (sqr(u) + sqr(v));

    numvector<numvector<double, dimPh>, dimPh> res;
    res =
    {
        // Rx
        {        1.0,  0.0,  0.0,   1.0,        1.0 },
        {     u - cS,  0.0,  0.0,     u,     u + cS },
        {          v,  1.0,  0.0,     v,          v },
        {        0.0,  0.0,  1.0,   0.0,        0.0 },
        { H - cS * u,    v,  0.0, magU2, H + cS * u }
    };

    return res;
}
numvector<numvector<double, dimPh>, dimPh> Physics::getRy(const numvector<double, dimPh>& sol) const
{
	double cS = c(sol);
	double rho = sol[0];
	double u = sol[1] / rho;
	double v = sol[2] / rho;
	double p = getPressure(sol);

	double H = (sol[4] + p) / rho;
	double magU2 = 0.5 * (sqr(u) + sqr(v));

	numvector<numvector<double, dimPh>, dimPh> res;
	res =
	{
		// Ry
		{ 1.0,			0.0,  0.0,   1.0,        1.0 },
		{ u,		   -1.0,  0.0,     u,          u },
		{ v - cS,		0.0,  0.0,     v,     v + cS },
		{ 0.0,			0.0,  1.0,   0.0,        0.0 },
		{ H - cS * v,   -u,   0.0, magU2, H + cS * v },
	};

	return res;
}

numvector<numvector<double, dimPh>, dimPh> Physics::getL(const numvector<double, dimPh>& sol, const Point& n) const
{
    double cS = c(sol);
    double rho = sol[0];
    double u = sol[1] / rho;
    double v = sol[2] / rho;

    double B1 = (cpcv - 1.0) / sqr(cS);
    double B2 = 0.5 * B1 * (sqr(u) + sqr(v));

    numvector<numvector<double, dimPh>, dimPh> res;

    res =
    {
        { 0.5 * (B2 + (u*n.x() + v*n.y())/cS), -0.5 * (B1 * u + n.x()/cS), -0.5 * (B1 * v + n.y()/cS),  0.0, 0.5 * B1},
        {                   u*n.y() - v*n.x(),                     -n.y(),                      n.x(),  0.0,      0.0},
        {                            1.0 - B2,                     B1 * u,                     B1 * v,  0.0,      -B1},
        {                                 0.0,                        0.0,                        0.0,  1.0,      0.0},
        { 0.5 * (B2 - (u*n.x() + v*n.y())/cS), -0.5 * (B1 * u - n.x()/cS), -0.5 * (B1 * v - n.y()/cS),  0.0, 0.5 * B1}
    };


    return res;
}
numvector<numvector<double, dimPh>, dimPh> Physics::getR(const numvector<double, dimPh>& sol, const Point& n) const
{
    double cS = c(sol);
    double rho = sol[0];
    double u = sol[1] / rho;
    double v = sol[2] / rho;
    double p = getPressure(sol);

    double H = (sol[4] + p) / rho;
    double magU2 = 0.5 * (sqr(u) + sqr(v));

    numvector<numvector<double, dimPh>, dimPh> res;

    res =
    {
        {                          1.0,                0.0,   1.0,  0.0,                          1.0 },
        {               u - cS * n.x(),             -n.y(),     u,  0.0,                 u + cS*n.x() },
        {               v - cS * n.y(),              n.x(),     v,  0.0,                 v + cS*n.y() },
        {                          0.0,                0.0,   0.0,  1.0,                          0.0 },
        { H - cS * (u*n.x() + v*n.y()), -u*n.y() + v*n.x(), magU2,  0.0, H + cS * (u*n.x() + v*n.y()) }
    };

    return res;
}









