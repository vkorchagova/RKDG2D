#include "Problem.h"
#include <iostream>

using namespace std;

// ------------------ Constructors & Destructors ----------------

Problem::Problem(const Mesh2D& mesh2D)
{
    mesh = mesh2D;

    //phi[0] = [](numvector<double,2> r){ return 0.5; };
    //phi[1] = [](numvector<double,2> r){ return 0.5*sqrt(3)*r[0]; };
    //phi[2] = [](numvector<double,2> r){ return 0.5*sqrt(3)*r[1]; };

    //gradPhi[0] = [](numvector<double,2> r){ return numvector<double,2>{ 0.0, 0.0 }; };
    //gradPhi[1] = [](numvector<double,2> r){ return numvector<double,2>{ 0.5*sqrt(3), 0.0 }; };
    //gradPhi[2] = [](numvector<double,2> r){ return numvector<double,2>{ 0.0, 0.5*sqrt(3) }; };
	
	double hx = mesh.hx;
	double hy = mesh.hy;
	
    phi[0] = [=](numvector<double, 2> r, int iCell){ return 1.0 / sqrt(hx*hy); };
    phi[1] = [=](numvector<double, 2> r, int iCell){ return sqrt(12.0 / hy / pow(hx, 3)) * (r[0] - mesh.cellCenters[iCell][0]); };
    phi[2] = [=](numvector<double, 2> r, int iCell){ return sqrt(12.0 / hx / pow(hy, 3)) * (r[1] - mesh.cellCenters[iCell][1]); };

	gradPhi[0] = [=](numvector<double, 2> r, int iCell){ return numvector<double, 2>{ 0.0, 0.0 }; };
	gradPhi[1] = [=](numvector<double, 2> r, int iCell){ return numvector<double, 2>{ sqrt(12.0 / hy / pow(hx, 3)), 0.0 }; };
	gradPhi[2] = [=](numvector<double, 2> r, int iCell){ return numvector<double, 2>{ 0.0, sqrt(12.0 / hx / pow(hy, 3)) }; };

    writer.open("alphaCoeffs");

} // end constructor by mesh

Problem::~Problem()
{
    writer.close();
}

numvector<double, 5> Problem::reconstructSolution(numvector<double, \
												nShapes + 5>& alpha, \
												numvector<double, 2>& point, \
												int iCell)
{
	numvector<double, 5> sol(0.0);

    for (int i = 0; i < 5; ++i)
	{
		for (int j = 0; j < nShapes; ++j)
			sol[i] += phi[j](point,iCell)*alpha[i*nShapes + j];
	}

	return sol;   
} // end reconstructSolution


// ------------------ Private class methods --------------------

void Problem::write(ostream& writer, const numvector<double,5*nShapes>& coeffs)
{
    for (int i = 0; i < 5*nShapes; ++i)
        writer << coeffs[i] << ' ';

    writer << endl;
}

double Problem::getPressure(numvector<double, 5> sol)
{
	double magU = pow(sol[2], 2) + pow(sol[3], 2) + pow(sol[1], 2);

	return (cpcv - 1)*(sol[4] - 0.5*magU / sol[0]);  
} // end getPressure

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

    double rho0 = 1;

    function<double(const numvector<double,2> r)> initRho = \
        [](const numvector<double,2> r) \
            { return 0.001 * exp(-2*pow(r[0]-2,2) - 2*pow(r[1]-2.5,2)); };

    function<double(const numvector<double,2>& r)> init[5];

    init[0] = [&](const numvector<double,2> r) { return rho0 + initRho(r); } ;
    init[1] = [](const numvector<double,2> r) { return 0.0; };
    init[2] = [](const numvector<double,2> r) { return 0.0; };
    init[3] = [](const numvector<double,2> r) { return 0.0; };
    init[4] = [&](const numvector<double,2> r) { return (rho0 + initRho(r)) / cpcv / (cpcv - 1.0); };


	int nCells = mesh.cells.size();

	alphaPrev.reserve(nCells);

    GaussIntegrator GP;

    // for internal cells
    for (int k = 0; k < mesh.nInternalCells; ++k)
	{
        numvector<numvector<double,2>,4> nodes = mesh.getCellCoordinates(k);

        for (int p = 0; p < 5; ++p)
		{
			for (int q = 0; q < nShapes; ++q)
			{
                function<double(const numvector<double,2>& x)> f = \
                        [=](const numvector<double,2>& x){return phi[q](x,k)*init[p](x);};

                alphaPrev[k][p * 3 + q] = GP.integrate(f, nodes);
			}
		}

       write(writer,alphaPrev[k]);
	}

    // for ghost cells
/*
    numvector<double, 5*nShapes> infCondition = {rho0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, (rho0) / cpcv / (cpcv - 1.0), 0.0, 0.0};

    for (int k = mesh.nInternalCells; k < nCells; ++k)
    {
        double sqrtJ = sqrt(mesh.hx * mesh.hy);

        alphaPrev[k] = sqrtJ * infCondition;

        write(writer,alphaPrev[k]);

    }
*/
} // end setInitialConditions


// ------------------ Public class methods --------------------
