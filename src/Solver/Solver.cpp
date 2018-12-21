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
    int nCells = M.cells.size();
    numvector<double, dimS> alpha;

    sln.SOL.resize(nCells);
    sln.SOL_aux.resize(nCells);

//#pragma omp parallel for
    for (int k = 0; k < nCells; ++k)
    {
        alpha = projection(prb.init, k);
        sln.SOL[k] = correctNonOrthoCell(alpha, k);
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

    int nCells = M.cells.size();

    numvector<double, dimS> rhs;

    for (int k = 0; k < nCells; ++k)
    {
        for (int j = 0; j < dimS; ++j)
            reader >> sln.SOL[k][j];
    }// for k

    reader.close();
}

numvector<double, dimS> Solver::projection(const std::function<numvector<double, dimPh>(const Point& point)>& foo, int iCell) const
{
	numvector<double, dimS> alpha;	
	numvector<double, dimPh> buffer;
	for (int q = 0; q < nShapes; ++q)
	{
		std::function<numvector<double, dimPh>(const Point&)> f = \
			[&](const Point& p) { return B.phi[q](iCell, p) * foo(p); };

		buffer = integrate(*(M.cells[iCell]), f);
		for (int p = 0; p < dimPh; ++p)
		{
			alpha[p*nShapes + q] = buffer[p];
		}// for p
	}// for shapes

	return alpha;
} // end projection


  /// REWRITE
vector<numvector<double, dimS>> Solver::assembleRHS(const std::vector<numvector<double, dimS>> &SOL)
{
    int nCells = M.cells.size();

    vector<numvector<double, dimS>> rhs(nCells);
  /*  double ts0 = 0;
    double te0 = 0;

    
    
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

    */

    return rhs;
}

numvector<double, dimS> Solver::correctNonOrthoCell(const numvector<double, dimS>& alpha, int iCell) const
{
	numvector<double, dimS> alphaCorr;

	vector<double> solution(nShapes); //for 3 ff!!!


	for (int iSol = 0; iSol < dimPh; ++iSol)
	{
		// solve slae
		solution[0] = alpha[iSol*nShapes];


		if (nShapes == 3)
		{
			solution[2] = (alpha[iSol*nShapes + 2] * B.gramian[iCell][0] - alpha[iSol*nShapes + 1] * B.gramian[iCell][1]) \
				            / (B.gramian[iCell][2] * B.gramian[iCell][0] - B.gramian[iCell][1] * B.gramian[iCell][1]);
			solution[1] = (alpha[iSol*nShapes + 1] - solution[2] * B.gramian[iCell][1]) \
				            / (B.gramian[iCell][0]);

            //cout << "\tiSol = " << iSol << endl;

            //solution[2] = (rhs[iSol*nShapes + 2] * gramian[1][1] - rhs[iSol*nShapes + 1] * gramian[2][1]) \
            //    / (gramian[2][2] * gramian[1][1] - gramian[1][2] * gramian[2][1]);
            //solution[1] = (rhs[iSol*nShapes + 1] - solution[2] * gramian[1][2]) \
            //      / (gramian[1][1]);
		}

		//set solution to appropriate positions
		for (int i = 0; i < nShapes; ++i)
        {
			//cout << "i = " << i << "; iAlpha = " << i + iSol*nShapes << endl; 
            alphaCorr[i + iSol*nShapes] = solution[i];
        }

        //cout << "the next sol..." << endl;
	}

    //cout << "result in fun = " << alphaCorr << endl;
    //cout << "the next cell..." << endl;

	return alphaCorr;
}


vector<numvector<double, dimS>> Solver::correctNonOrtho(const vector<numvector<double, dimS>>& alpha) const
{
	vector<numvector<double, dimS>> alphaCorr(alpha.size());

	cout << "size of alpha " << alpha.size() << endl;
    cout << "size of alphaCorr " << alphaCorr.size() << endl;
    for (size_t iCell = 0; iCell < M.nRealCells; ++iCell)
	{
		cout << "iCell in common  = " << iCell << endl;
        alphaCorr[iCell] = correctNonOrthoCell(alpha[iCell], iCell);
        cout << "result = " << alphaCorr[iCell] << endl;
	}// for cells

	return alphaCorr;
}


numvector<double, dimS> Solver::correctPrevIterCell(const numvector<double, dimS>& alphaCorr, int iCell) const
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


vector<numvector<double, dimS>> Solver::correctPrevIter() const
{
	vector<numvector<double, dimS>> alpha;
	alpha.resize(sln.SOL.size());

	//for(int cell=0;)
	//   alpha[cell]... SOL OR SOL_aux

	return alpha;
}