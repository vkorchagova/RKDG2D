#include "Problem.h"


Problem::Problem()
{
}

Problem::Problem(Mesh2D& mesh2D)
{
	mesh = mesh2D;

	phi.reserve(nShapes);

	phi.push_back(phi0);
	phi.push_back(phi1);
	phi.push_back(phi2);
}


Problem::~Problem()
{
}

numvector<double, 5> Problem::reconstructSolution(numvector<double, \
												nShapes + 5>& alpha, \
												numvector<double, 2>& localCoord)
{
	numvector<double, 5> sol(0.0);

	for (int i = 0; i < 5; ++i)
	{
		for (int j = 0; j < nShapes; ++j)
			sol[i] += phi[j](localCoord)*alpha[i+5*j];
	}

	return sol;
}

double Problem::getPressure(numvector<double, 5> sol)
{
	double magU = pow(sol[2], 2) + pow(sol[3], 2) + pow(sol[1], 2);

	return (gamma - 1)*(sol[4] - 0.5*magU / sol[0]);
}

numvector<double, 5> Problem::fluxF(numvector<double, 5> sol)
{
	double u = sol[1] / sol[0];
	double p = getPressure(sol);

	return{ sol[1], u*sol[1] + p, u*sol[2], u*sol[3], (sol[4] + p)*u };
}

numvector<double, 5> Problem::fluxG(numvector<double, 5> sol)
{
	double v = sol[2] / sol[0];
	double p = getPressure(sol);

	return{ sol[2], v*sol[1], v*sol[2] + p, v*sol[3], (sol[4] + p)*v };
}
