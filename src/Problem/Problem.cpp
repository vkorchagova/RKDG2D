#include "Problem.h"
#include <iostream>

using namespace std;

// ------------------ Constructors & Destructors ----------------

Problem::Problem(Mesh2D &mesh2D)
{
    mesh = &mesh2D;

    writer.open("alphaCoeffs");

} // end constructor by mesh

Problem::~Problem()
{
    writer.close();
}


// ------------------ Private class methods --------------------

// ------------------ Public class methods --------------------

void Problem::write(ostream& writer, const numvector<double,5*nShapes>& coeffs)
{
    for (int i = 0; i < 5*nShapes; ++i)
        writer << coeffs[i] << ' ';

    writer << endl;
} //end write

double Problem::getPressure(numvector<double, 5> sol)
{
	double magU = pow(sol[2], 2) + pow(sol[3], 2) + pow(sol[1], 2);

	return (cpcv - 1)*(sol[4] - 0.5*magU / sol[0]);  
} // end getPressure

double Problem::c(numvector<double, 5> sol)
{
    return sqrt( cpcv * getPressure(sol) / sol[0]);
} // end c for cell

double Problem::c_av(numvector<double, 5> solOne, numvector<double, 5> solTwo)
{
    double semiRho = 0.5*(solOne[0] + solTwo[0]);
    double semiP =   0.5*(getPressure(solOne) + getPressure(solTwo));

    return sqrt( cpcv * semiP / semiRho);
} // end c for edge

numvector<double, 5> Problem::lambdaF(numvector<double, 5> solOne, numvector<double, 5> solTwo)
{
    double u = fabs(0.5*(solOne[1] / solOne[0] + solTwo[1] / solTwo[0])); //estimation!!!
    double soundSpeed = c_av(solOne, solTwo);

    return {u - soundSpeed, u, u, u, u + soundSpeed};
} // end lambdaF

numvector<double, 5> Problem::lambdaG(numvector<double, 5> solOne, numvector<double, 5> solTwo)
{
    double v = fabs(0.5*(solOne[2] / solOne[0] + solTwo[2] / solTwo[0])); //estimation!!!
    double soundSpeed = c_av(solOne, solTwo);

    return {v - soundSpeed, v, v, v, v + soundSpeed};
} // end lambdaG


numvector<double, 5> Problem::fluxF(numvector<double, 5> sol)
{
	double u = sol[1] / sol[0];
	double p = getPressure(sol);

    return { sol[1], u*sol[1] + p, u*sol[2], u*sol[3], (sol[4] + p)*u };
} // end fluxF

numvector<double, 5> Problem::fluxG(numvector<double, 5> sol)
{
	double v = sol[2] / sol[0];
	double p = getPressure(sol);

    return { sol[2], v*sol[1], v*sol[2] + p, v*sol[3], (sol[4] + p)*v };
} // end fluxG

void Problem::setInitialConditions()
{
    // define functions for initial conditions

    double rho0 = 1.0;

    function<double(const Point& r)> initRho = \
        [](const Point& r) \
            { return 0.001 * exp( -2.0 * pow(r.x() - 2.0, 2) - 2.0 * pow(r.y() - 2.0, 2)); };

    function<numvector<double, 5>(const Point& r)> init = \
        [&](const Point& r) { return numvector<double, 5> { rho0 + initRho(r), 0.0, 0.0, 0.0, (rho0 + initRho(r)) / cpcv / (cpcv - 1.0) }; };


    int nCells = mesh->nInternalCells;

    //alphaPrev.reserve(nCells);

    //GaussIntegrator GP;
    //numvector<double, 5> buffer;

    // for internal cells
    for (int k = 0; k < nCells; ++k)
	{
        mesh->cells[k].setProblem(*this);
        mesh->cells[k].setLocalInitialConditions(init);

        //write(writer,mesh->cells[k].alphaPrev);
	}
    // for ghost cells
/*
    numvector<double, 5*nShapes> infCondition = {rho0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, (rho0) / cpcv / (cpcv - 1.0), 0.0, 0.0};

    for (int k = mesh->nInternalCells; k < nCells; ++k)
    {
        double sqrtJ = sqrt(mesh->hx * mesh->hy);

        alphaPrev[k] = sqrtJ * infCondition;

        //write(writer,alphaPrev[k]);

    }
*/

} // end setInitialConditions







