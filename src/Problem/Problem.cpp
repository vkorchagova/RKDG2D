#include "Problem.h"
#include <iostream>

using namespace std;

// ------------------ Constructors & Destructors ----------------

Problem::Problem()
{
    // Function for initial condition
//    double rho0 = 1.0;
//    double e0 = rho0  / cpcv / (cpcv - 1.0) ;

//    function<double(const Point& r)> initRho = [](const Point& r) \
//            { return 0.001 * exp( -2.0 * pow(r.x() - 2.0, 2) - 2.0 * pow(r.y() - 2.0, 2)); };

//    init = [=](const Point& r) { return numvector<double, 5> { rho0 + initRho(r), 0.0, 0.0, 0.0, (rho0 + initRho(r)) / cpcv / (cpcv - 1.0) }; };


    // For boundary conditions
//    infty = {rho0, 0.0, 0.0, 0.0, e0};

} // end constructor by mesh

Problem::~Problem()
{

}


// ------------------ Private class methods --------------------

// ------------------ Public class methods --------------------



//// RKDG methods

void Problem::getAlpha(const std::vector<numvector<double, 5 * nShapes> >& a)
{
    alpha = a;
}

double Problem::getPressure(const numvector<double, 5>& sol) const
{
	double magU = pow(sol[2], 2) + pow(sol[3], 2) + pow(sol[1], 2);

	return (cpcv - 1)*(sol[4] - 0.5*magU / sol[0]);
} // end getPressure

double Problem::c(const numvector<double, 5>& sol) const
{
    return sqrt( cpcv * getPressure(sol) / sol[0]);
} // end c for cell

double Problem::c_av(const numvector<double, 5>& solOne, const numvector<double, 5>& solTwo) const
{
    double semiRho = 0.5*(solOne[0] + solTwo[0]);
    double semiP =   0.5*(getPressure(solOne) + getPressure(solTwo));

    return sqrt( cpcv * semiP / semiRho);
} // end c for edge

numvector<double, 5> Problem::lambdaF(const numvector<double, 5>& solOne, const numvector<double, 5>& solTwo) const
{
    double u = fabs(0.5*(solOne[1] / solOne[0] + solTwo[1] / solTwo[0])); //estimation!!!
    double soundSpeed = c_av(solOne, solTwo);

    return {u - soundSpeed, u, u, u, u + soundSpeed};
} // end lambdaF

numvector<double, 5> Problem::lambdaG(const numvector<double, 5>& solOne, const numvector<double, 5>& solTwo) const
{
    double v = fabs(0.5*(solOne[2] / solOne[0] + solTwo[2] / solTwo[0])); //estimation!!!
    double soundSpeed = c_av(solOne, solTwo);

    return {v - soundSpeed, v, v, v, v + soundSpeed};
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









