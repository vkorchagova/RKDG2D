#include "Solution.h"

#include <iostream>

using namespace std;


numvector<double, dimPh> Solution::reconstructSolution(const Point& point ) const
{
    numvector<double, dimPh> sol(0.0);
    vector<double> phip(nShapes);

#pragma omp simd
    for (int j = 0; j < nShapes; ++j)
        phip[j] = phi[j](point);

    for (int i = 0; i < dimPh; ++i)
    {
#pragma omp simd
        for (int j = 0; j < nShapes; ++j)
            sol[i] += phip[j] * SOL[number][i * nShapes + j];
    }

    return sol;

} // end reconstructSolution
numvector<double, dimPh> Solution::reconstructSolution(const Point& point, const numvector<double, dimS>& SOL) const
{
	numvector<double, dimPh> sol(0.0);
	vector<double> phip(nShapes);

#pragma omp simd
	for (int j = 0; j < nShapes; ++j)
		phip[j] = phi[j](point);

	for (int i = 0; i < dimPh; ++i)
	{

#pragma omp simd
		for (int j = 0; j < nShapes; ++j)
			sol[i] += phip[j] * SOL[i * nShapes + j];
	}

	return sol;

} // end reconstructSolution

double Solution::reconstructSolution(const Point& point, Variables var ) const
{
    double sol(0.0);

#pragma omp simd    
    for (int j = 0; j < nShapes; ++j)
        sol += phi[j](point) * SOL[number][var * nShapes + j];

    return sol;
} // end reconstructSolution
double Solution::reconstructSolution(const Point& point, const numvector<double, dimS>& SOL, Variables var ) const
{
    double sol(0.0);

#pragma omp simd
    for (int j = 0; j < nShapes; ++j)
        sol += phi[j](point) * SOL[var * nShapes + j];

    return sol;
} // end reconstructSolution

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
void Solution::writeSolutionVTK(string fileName) const
{
	ofstream output;
	output.precision(16);
	output.open(fileName + ".vtk");

	mesh.exportMeshVTK_polyvertices(output);

	// get cell data

	output << "CELL_DATA " << mesh.nCells << endl;

	output << "SCALARS rho double" << endl;
	output << "LOOKUP_TABLE default" << endl;

	for (const shared_ptr<Cell> cell : mesh.cells)
		output << cell->reconstructSolution(cell->getCellCenter(), 0) << endl;

	output << "SCALARS e double" << endl;
	output << "LOOKUP_TABLE default" << endl;

	for (const shared_ptr<Cell> cell : mesh.cells)
		output << cell->reconstructSolution(cell->getCellCenter(), 4) << endl;

	output << "SCALARS p double" << endl;
	output << "LOOKUP_TABLE default" << endl;

	for (const shared_ptr<Cell> cell : mesh.cells)
		output << problem.getPressure(cell->reconstructSolution(cell->getCellCenter())) << endl;


	output << "VECTORS U double" << endl;

	for (const shared_ptr<Cell> cell : mesh.cells)
	{
		double rho = cell->reconstructSolution(cell->getCellCenter(), 0);
		output << cell->reconstructSolution(cell->getCellCenter(), 1) / rho << endl;
		output << cell->reconstructSolution(cell->getCellCenter(), 2) / rho << endl;
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
	output << cell->reconstructSolution(cell->nodes[j], 0) << endl;

	output << "SCALARS e double" << endl;
	output << "LOOKUP_TABLE default" << endl;

	for (const shared_ptr<Cell> cell : mesh.cells)
	for (int j = 0; j < cell->nEntities; ++j)
	output << cell->reconstructSolution(cell->nodes[j], 4) << endl;

	output << "SCALARS p double" << endl;
	output << "LOOKUP_TABLE default" << endl;

	for (const shared_ptr<Cell> cell : mesh.cells)
	for (int j = 0; j < cell->nEntities; ++j)
	output << problem.getPressure(cell->reconstructSolution(cell->nodes[j])) << endl;


	output << "VECTORS U double" << endl;

	for (const shared_ptr<Cell> cell : mesh.cells)
	for (int j = 0; j < cell->nEntities; ++j)
	{
	double rho = cell->reconstructSolution(cell->nodes[j], 0);
	output << cell->reconstructSolution(cell->nodes[j], 1) / rho << endl;
	output << cell->reconstructSolution(cell->nodes[j], 2) / rho << endl;
	output << 0.0 << endl;
	}

	*/
	output.close();
}


