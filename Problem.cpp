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

} // end constructor by mesh

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
            { return 0.001 * exp(-2*pow(r[0] - 0.5,2) - 2*pow(r[1]-0.5,2)); };

    function<double(const numvector<double,2>& r)> init[5];

    init[0] = [&](const numvector<double,2> r) { return rho0 + initRho(r); } ;
    init[1] = [](const numvector<double,2> r) { return 0.0; };
    init[2] = [](const numvector<double,2> r) { return 0.0; };
    init[3] = [](const numvector<double,2> r) { return 0.0; };
    init[4] = [&](const numvector<double,2> r) { return (rho0 + initRho(r)) / cpcv / (cpcv - 1.0); };

    // get start alpha coeffs
    //numvector < numvector<double, 2>,4> gPointsGlobal;
    //numvector < double, nGPoints> fValues;

	

    //double J = mesh.hx * mesh.hy * 0.25;

	int nCells = mesh.cells.size();

	alphaPrev.reserve(nCells);

    GaussIntegrator GP;

    for (int k = 0; k < mesh.nInternalCells; ++k)
	{
        numvector<numvector<double,2>,4> nodes = mesh.getCellCoordinates(k);
        cout << nodes << endl;

        for (int p = 0; p < 5; ++p)
		{
			for (int q = 0; q < nShapes; ++q)
			{
                //for (int i = 0; i < nGPoints; ++i)
                //{
                //	gPointsGlobal = mesh.localToGlobal(k, gPoints[i]);
                //	fValues[i] = phi[q](gPointsGlobal[i],k) * init[p](gPointsGlobal[i]);

					//cout << "fVal = " << fValues[i] << ' ';
                //}

				//cout << endl;

                function<double(const numvector<double,2>& x)> f = \
                        [=](const numvector<double,2>& x){return phi[q](x,k)*init[p](x);};

                //alphaPrev[k][p * 3 + q] = integrate2D(fValues, J);
                alphaPrev[k][p * 3 + q] = GP.integrate(f, nodes);
			}
		}
	}

} // end setInitialConditions


// ------------------ Public class methods --------------------
