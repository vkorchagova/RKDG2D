#include <iostream>
#include <omp.h>

#include "Physics.h"
#include "Params.h"

using namespace std;

// ------------------ Constructors & Destructors ----------------

Physics::Physics()
{
    cpcv = 1.4; // default value
    covolume = 0.0; // default value
    R = 1.0/cpcv;
} // end constructor 

Physics::~Physics()
{

}


// ------------------ Private class methods --------------------

// ------------------ Public class methods --------------------


double Physics::magU(const numvector<double, dimPh>& sol) const
{
    double rhoV2 = sqr(sol[1]) + sqr(sol[2]) + sqr(sol[3]);

    return sqrt(rhoV2) / sol[0];
} // end magU


double Physics::getPressure(const numvector<double, dimPh>& sol) const
{
    double rhoV2 = sqr(sol[1]) + sqr(sol[2]) + sqr(sol[3]);
    //double rhoEps = sol[4] - 0.2 * rhoV2 / sol[0];

    //double p = (cpcv - 1.0)*(sol[4] - 0.5*rhoV2 / sol[0]);  // for ideal gas

    double p = (cpcv - 1.0)*(sol[4] - 0.5*rhoV2 / sol[0]) / (1.0 - sol[0] * covolume); 

    return p;
} // end getPressure


double Physics::getMachNumber(const numvector<double, dimPh>& sol) const
{
    double rhoV2 = sqr(sol[1]) + sqr(sol[2]) + sqr(sol[3]);

    return sqrt(rhoV2) / sol[0] / this->c(sol);
} // end getPressure


double Physics::c(const numvector<double, dimPh>& sol) const
{
    //double c2 = cpcv * getPressure(sol) / sol[0];

    double c2 = cpcv * getPressure(sol) / sol[0] / (1.0 - sol[0] * covolume);

    if (c2 < 0)
    {
        cout << "sound speed < 0 =( !!!! State = " << sol << endl;
        exit(1);
    }

    return sqrt( c2 );
} // end c for cell


double Physics::e(double rho, double u, double v, double w, double p) const
{
    return p / (cpcv - 1.0) * (1.0 - rho * covolume) + 0.5 * rho * (u*u + v*v + w*w);
} // end e for cell


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


void Physics::computeLR(
        const numvector<double, dimPh>& sol, 
        const Point& n, 
        numvector<numvector<double, dimPh>, dimPh>& L,
        numvector<numvector<double, dimPh>, dimPh>& R
    ) const
{
    double cS = c(sol);
    double rho = sol[0];
    double u = sol[1] / rho;
    double v = sol[2] / rho;
    double magU2 = 0.5 * (sqr(u) + sqr(v));

    double B1 = (cpcv - 1.0) / sqr(cS);
    double B2 = B1 * magU2;

    double p = getPressure(sol);

    double H = (sol[4] + p) / rho;

    double Un = u*n.x() + v*n.y();
    double Unort = u*n.y() - v*n.x();
    double UndcS = Un/cS;
    double nxdcS = n.x()/cS;
    double nydcS = n.y()/cS;
    double B1u = B1 * u;
    double B1v = B1 * v;
    double cSnx = cS * n.x();
    double cSny = cS * n.y();
    double cSUn = cS * Un;

    L =
    {
        { 0.5 * (B2 + UndcS), -0.5 * (B1u + nxdcS), -0.5 * (B1v + nydcS),  0.0, 0.5 * B1},
        {              Unort,               -n.y(),                n.x(),  0.0,      0.0},
        {           1.0 - B2,                  B1u,                  B1v,  0.0,      -B1},
        {                0.0,                  0.0,                  0.0,  1.0,      0.0},
        { 0.5 * (B2 - UndcS), -0.5 * (B1u - nxdcS), -0.5 * (B1v - nydcS),  0.0, 0.5 * B1}
    };

    R =
    {
        {      1.0,    0.0,   1.0,  0.0,      1.0 },
        { u - cSnx, -n.y(),     u,  0.0, u + cSnx },
        { v - cSny,  n.x(),     v,  0.0, v + cSny },
        {      0.0,    0.0,   0.0,  1.0,      0.0 },
        { H - cSUn, -Unort, magU2,  0.0, H + cSUn }
    };
}










