#include "Problem.h"
#include <iostream>

using namespace std;

// ------------------ Constructors & Destructors ----------------

Problem::Problem()
{
    // Function for initial condition
    double rho0 = 1.0;
    double e0 = rho0  / cpcv / (cpcv - 1.0) ;
    double v0 = 0.0;

    function<double(const Point& r)> initRho = [=](const Point& r) \
    { 
    //    return 1.0;
    //    return rho0 + 0.001 * exp( -2.0 * sqr(r.x() - 2.0) - 2.0 * sqr(r.y() - 2.0));
    //    return (r.y() < 0.5) ? 1.0 : 0.125;
        return (r.x() < 0.5) ? 1.0 : 0.125;
    //   return (r.y() < 1.0 && r.x() < 1.0 && r.y() > 2.0 && r.x() > 2.0) ? 0.0 : 1.0;
    //return (r.y() < 0.5) ? r.y() + 0.01 : r.y() + 0.51;
    };

    function<double(const Point& r)> initP = [=](const Point& r) \
    { 
    //    return (initRho(r)) / cpcv;
    //    return (r.y() < 0.5) ? 1.0 : 0.1;
        return (r.x() < 0.5) ? 1.0 : 0.1;
    };


    function<double(const Point& r)> initV = [](const Point& r) \
    {
    //    return 0.0;
    //    return (r.y() < 0.5) ? 0.0 : 0.0;
        return (r.x() < 0.5) ? 0.0 : 0.0;
	//return 0.001 * exp( -2.0 * pow(r.x() - 4.0, 2) - 2.0 * pow(r.y() - 4.0, 2)); 
    //return r.y();
    };

    //init = [=](const Point& r) { return numvector<double, 5> { initRho(r), 0.0, initV(r), 0.0, initP(r) / (cpcv - 1.0) }; };
    init = [=](const Point& r) { return numvector<double, 5> { initRho(r), initV(r), 0.0, 0.0, initP(r) / (cpcv - 1.0) }; };


    // For boundary conditions
    infty = { rho0, 0.0, v0, 0.0, e0 };

} // end constructor by mesh

Problem::~Problem()
{

}


// ------------------ Private class methods --------------------

// ------------------ Public class methods --------------------



//// RKDG methods

void Problem::setAlpha(const std::vector<numvector<double, 5 * nShapes> >& a)
{
    alpha = a;
}

double Problem::getPressure(const numvector<double, 5>& sol) const
{
	double magRhoU2 = sqr(sol[1]) + sqr(sol[2]) + sqr(sol[3]);

	return (cpcv - 1.0)*(sol[4] - 0.5*magRhoU2 / sol[0]);
} // end getPressure

double Problem::c(const numvector<double, 5>& sol) const
{
    return sqrt( cpcv * getPressure(sol) / sol[0]);
} // end c for cell

//double Problem::h(const numvector<double, 5>& sol) const
//{
//    return (sol[4] + getPressure(sol)) / sol[0];
//}

double Problem::c_av(const numvector<double, 5>& solOne, const numvector<double, 5>& solTwo) const
{
    double semiRho = 0.5*(solOne[0] + solTwo[0]);
    double semiP =   0.5*(getPressure(solOne) + getPressure(solTwo));

    return sqrt( cpcv * semiP / semiRho);
} // end c for edge

numvector<double, 5> Problem::lambdaF_Roe(const numvector<double, 5>& solOne, const numvector<double, 5>& solTwo) const
{
    double sqrtRhoLeft = sqrt(solOne[0]);
    double sqrtRhoRight = sqrt(solTwo[0]);
    double sumSqrtRho = sqrtRhoLeft + sqrtRhoRight;

    double u_av = ( solOne[1] / sqrtRhoLeft + solTwo[1] / sqrtRhoRight ) / sumSqrtRho;
    double h_av = ( (solOne[4] + getPressure(solOne)) / sqrtRhoLeft + (solTwo[4] + getPressure(solTwo)) / sqrtRhoRight ) / sumSqrtRho;

    double c_av = sqrt( (cpcv - 1) * (h_av - 0.5 * sqr(u_av)) );

    return {u_av - c_av, u_av, u_av, u_av, u_av + c_av};

}

numvector<double, 5> Problem::lambdaF_semisum(const numvector<double, 5>& solOne, const numvector<double, 5>& solTwo) const
{
    double u = 0.5*(solOne[1] / solOne[0] + solTwo[1] / solTwo[0]);
    double soundSpeed = c_av(solOne, solTwo);

    return { u - soundSpeed, u, u, u, u + soundSpeed};
} // end lambdaF

numvector<double, 5> Problem::lambdaF(const numvector<double, 5>& solOne, const numvector<double, 5>& solTwo) const
{
    return lambdaF_Roe(solOne,solTwo);
} // end lambdaF

numvector<double, 5> Problem::lambdaG(const numvector<double, 5>& solOne, const numvector<double, 5>& solTwo) const
{
    double v = 0.5*(solOne[2] / solOne[0] + solTwo[2] / solTwo[0]);
    double soundSpeed = c_av(solOne, solTwo);

    return { v - soundSpeed, v, v, v, v + soundSpeed};
} // end lambdaG


numvector<double, 5> Problem::fluxF(const numvector<double, 5>& sol) const
{
    double u = sol[1] / sol[0];
    double p = getPressure(sol);

    return { sol[1], u*sol[1] + p, u*sol[2], u*sol[3], (sol[4] + p)*u };
} // end fluxF

numvector<double, 5> Problem::fluxG(const numvector<double, 5>& sol) const
{
    double v = sol[2] / sol[0];
    double p = getPressure(sol);

    return { sol[2], v*sol[1], v*sol[2] + p, v*sol[3], (sol[4] + p)*v };
} // end fluxG









