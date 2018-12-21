#include "Solution.h"

#include <iostream>

using namespace std;

///

Solution::Solution(Basis& Bas) : B(Bas)
{

}

///

numvector<double, dimPh> Solution::reconstruct(int iCell, const Point& point ) const
{
    numvector<double, dimPh> sol(0.0);
    vector<double> phip(nShapes);

#pragma omp simd
    for (int j = 0; j < nShapes; ++j)
        phip[j] = B.phi[j](iCell, point);

    for (int i = 0; i < dimPh; ++i)
    {
#pragma omp simd
        for (int j = 0; j < nShapes; ++j)
            sol[i] += phip[j] * SOL[iCell][i * nShapes + j];
    }

    return sol;

} // end reconstruct
numvector<double, dimPh> Solution::reconstruct(int iCell, const Point& point, const numvector<double, dimS>& SOL) const
{
	numvector<double, dimPh> sol(0.0);
	vector<double> phip(nShapes);

#pragma omp simd
	for (int j = 0; j < nShapes; ++j)
		phip[j] = B.phi[j](iCell, point);

	for (int i = 0; i < dimPh; ++i)
	{

#pragma omp simd
		for (int j = 0; j < nShapes; ++j)
			sol[i] += phip[j] * SOL[i * nShapes + j];
	}

	return sol;

} // end reconstruct

double Solution::reconstruct(int iCell, const Point& point, Variables var ) const
{
    double sol(0.0);

#pragma omp simd    
    for (int j = 0; j < nShapes; ++j)
        sol += B.phi[j](iCell, point) * SOL[iCell][var * nShapes + j];

    return sol;
} // end reconstruct
double Solution::reconstruct(int iCell, const Point& point, const numvector<double, dimS>& SOL, Variables var ) const
{
    double sol(0.0);

#pragma omp simd
    for (int j = 0; j < nShapes; ++j)
        sol += B.phi[j](iCell, point) * SOL[var * nShapes + j];

    return sol;
} // end reconstruct

/*
void Solution::write(string fileName, const vector<numvector<double, dimS>>& coeffs) const
{
	ofstream writer;
	writer.open(fileName);
	writer.precision(16);

	for (size_t k = 0; k < coeffs.size(); ++k)
	{
		for (int i = 0; i < dimS; ++i)
			writer << coeffs[k][i] << ' ';

		writer << endl;
	}
	writer.close();

} //end write


void Solution::writeSolutionVTK(const Mesh& mesh, string fileName) const
{
	ofstream output;
	output.precision(16);
	output.open(fileName + ".vtk");

	exportMeshVTK_polyvertices(output);

	// get cell data

	output << "CELL_DATA " << mesh.nRealCells << endl;

	output << "SCALARS rho double" << endl;
	output << "LOOKUP_TABLE default" << endl;

//	for (const shared_ptr<Cell> cell : mesh.cells)
	for (int iCell = 0; iCell < mesh.nRealCells; ++iCell)
		output << reconstruct(iCell, mesh.cells[iCell]->getCellCenter(), 0) << endl;

	output << "SCALARS e double" << endl;
	output << "LOOKUP_TABLE default" << endl;

	for (int iCell = 0; iCell < mesh.nRealCells; ++iCell)
		output << reconstruct(iCell, mesh.cells[iCell]->getCellCenter(), 4) << endl;

	output << "SCALARS p double" << endl;
	output << "LOOKUP_TABLE default" << endl;

	for (int iCell = 0; iCell < mesh.nRealCells; ++iCell)
		output << problem.getPressure(cell->reconstruct(cell->getCellCenter())) << endl;


	output << "VECTORS U double" << endl;

	for (const shared_ptr<Cell> cell : mesh.cells)
	{
		double rho = cell->reconstruct(cell->getCellCenter(), 0);
		output << cell->reconstruct(cell->getCellCenter(), 1) / rho << endl;
		output << cell->reconstruct(cell->getCellCenter(), 2) / rho << endl;
		output << 0.0 << endl;
	}

	//indicator.writeTroubledCellsVTK();

	// get point data

	output << "POINT_DATA " << mesh.nEntitiesTotal << endl;
	/*
	output << "SCALARS rho double" << endl;
	output << "LOOKUP_TABLE default" << endl;

	for (const shared_ptr<Cell> cell : mesh.cells)
	for (int j = 0; j < cell->nEntities; ++j)
	output << cell->reconstruct(cell->nodes[j], 0) << endl;

	output << "SCALARS e double" << endl;
	output << "LOOKUP_TABLE default" << endl;

	for (const shared_ptr<Cell> cell : mesh.cells)
	for (int j = 0; j < cell->nEntities; ++j)
	output << cell->reconstruct(cell->nodes[j], 4) << endl;

	output << "SCALARS p double" << endl;
	output << "LOOKUP_TABLE default" << endl;

	for (const shared_ptr<Cell> cell : mesh.cells)
	for (int j = 0; j < cell->nEntities; ++j)
	output << problem.getPressure(cell->reconstruct(cell->nodes[j])) << endl;


	output << "VECTORS U double" << endl;

	for (const shared_ptr<Cell> cell : mesh.cells)
	for (int j = 0; j < cell->nEntities; ++j)
	{
	double rho = cell->reconstruct(cell->nodes[j], 0);
	output << cell->reconstruct(cell->nodes[j], 1) / rho << endl;
	output << cell->reconstruct(cell->nodes[j], 2) / rho << endl;
	output << 0.0 << endl;
	}

	
	output.close();
}

*/

