#include "Solver.h"
#include <iostream>
#include <omp.h>

using namespace std;

Solver::Solver( Basis& Bas, Mesh& msh, Solution& soln,
				Problem &prob, Flux& flx) : B(Bas), M(msh), sln(soln),
											prb(prob), flux(flx)
{
    
}

void Solver::setInitialConditions()
{
    int nCells = M.nCells;
    numvector<double, dimS> alpha;

//#pragma omp parallel for
    for (int k = 0; k < nCells; ++k)
    {
        // REWRITE!!!
		std::function<numvector<double, dimPh>(const Point& r)> init;
		alpha = projection(prb.init, M.cells[k]);
        sln.SOL[k] = correctNonOrtho(alpha, M.cells[k]);
    }
    
    cout << "OK"<<endl;

} // end setInitialConditions

void Solver::setDefinedCoefficients(string fileName)
{
    ifstream reader;
    reader.open(fileName);

    if (!reader.is_open())
    {
        cout << "File " << fileName << " is not found\n";
        exit(0);
    }// if open

    int nCells = M.nCells;

    numvector<double, dimS> rhs;

    for (int k = 0; k < nCells; ++k)
    {
        for (int j = 0; j < dimS; ++j)
            reader >> sln.SOL[k][j];
    }// for k

    reader.close();
}

numvector<double, dimS> Solver::projection(const std::function<numvector<double, dimPh>(const Point& point)>& foo, const Cell& cell) const
{
	numvector<double, dimS> alpha;	
	numvector<double, dimPh> buffer;
	for (int q = 0; q < nShapes; ++q)
	{
		std::function<numvector<double, dimPh>(const Point&)> f = \
			[&](const Point& p) {  return B.phi[q](cell, p) * foo(p); };

		buffer = integrate(cell, f);
		for (int p = 0; p < dimPh; ++p)
		{
			alpha[p*nShapes + q] = buffer[p];
		}// for p
	}// for shapes

	return alpha;
} // end projection

vector<numvector<double, dimS>> Solver::assembleRHS(const std::vector<numvector<double, dimS>> &SOL)
{
    double ts0 = 0;
    double te0 = 0;

    int nCells = M.nCells;
    
    ts0 = omp_get_wtime();
    
    // compute fluxes in gauss points on edges
#pragma omp parallel for \
shared (mesh) \
default(none)
    for (size_t i = 0; i < M.nEdges; ++i)
        /mesh.edges[i]->getLocalFluxes(flux);

    te0 = omp_get_wtime();
    
    //cout << "get fluxes " << te0 - ts0 << endl;

    ts0 = omp_get_wtime();

    vector<numvector<double, dimS>> rhs(nCells);

#pragma omp parallel for \
shared (nCells, rhs, mesh) \
default(none)
    for (int k = 0; k < nCells; ++k) // for all cells
    {
        // compute internal integral
        /rhs[k] = mesh.cells[k]->cellIntegral();

        // compute boundary integrals
#pragma omp simd
        for (int i = 0; i < /mesh.cells[k]->nEntities; ++i)
            /rhs[k] -= mesh.cells[k]->edges[i]->boundaryIntegral(mesh.cells[k]);
    }// for cells

    te0 = omp_get_wtime();
    
    //cout << "get rhs " << te0 - ts0 << endl;

    return rhs;
}

numvector<double, dimS> Solver::correctNonOrtho(const numvector<double, dimS>& alpha, const Cell& cell) const
{
	numvector<double, dimS> alphaCorr;

	vector<double> solution(nShapes); //for 3 ff!!!

	for (int iSol = 0; iSol < dimPh; ++iSol)
	{
		// solve slae
		solution[0] = alpha[iSol*nShapes] / B.gramian[0][0];

		if (nShapes == 3)
		{
			solution[2] = (alpha[iSol*nShapes + 2] * B.gramian[1][1] - alpha[iSol*nShapes + 1] * B.gramian[2][1]) \
				/ (B.gramian[2][2] * B.gramian[1][1] - B.gramian[1][2] * B.gramian[2][1]);
			solution[1] = (alpha[iSol*nShapes + 1] - solution[2] * B.gramian[1][2]) \
				/ (B.gramian[1][1]);
		}

		//set solution to appropriate positions
		for (int i = 0; i < nShapes; ++i)
			alphaCorr[i + iSol*nShapes] = solution[i];
	}

	return alphaCorr;
}

numvector<double, dimS> Solver::correctPrevIter(const numvector<double, dimS>& alphaCorr, const Cell& cell) const
{
	numvector<double, dimS> alpha(0.0);
	for (int iSol = 0; iSol < dimPh; ++iSol)
	{
		for (int i = 0; i < nShapes; ++i)
#pragma omp simd
			for (int j = 0; j < nShapes; ++j)
				alpha[i + iSol*nShapes] += B.gramian[i][j] * alphaCorr[iSol*nShapes + j];
	}// for variables

	return alpha;
}