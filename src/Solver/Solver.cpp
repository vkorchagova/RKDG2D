#include "Solver.h"
#include <iostream>

using namespace std;

Solver::Solver( Mesh2D& msh, Problem &prb, Flux& flx) : mesh(msh),problem(prb),flux(flx)
{
    alphaPrev.resize(mesh.nCells);
    alphaNext.resize(mesh.nCells);
}


void Solver::write(string fileName, const vector<numvector<double,5*nShapes>>& coeffs) const
{
    ofstream writer;
    writer.open(fileName);
    writer.precision(16);

    for (size_t k = 0; k < coeffs.size(); ++k)
    {
        for (int i = 0; i < 5*nShapes; ++i)
            writer << coeffs[k][i] << ' ';

        writer << endl;
    }
    writer.close();

} //end write

void Solver::writeSolutionVTK(string fileName) const
{
    ofstream output;
    output.precision(16);
    output.open(fileName + ".vtk");

    mesh.exportMeshVTK_polyvertices(output);

    // get cell data

    output << "CELL_DATA " << mesh.nCells << endl;

//    output << "SCALARS rho double" << endl;
//    output << "LOOKUP_TABLE default" << endl;

//    for (const shared_ptr<Cell> cell : mesh.cells)
//        output << cell->reconstructSolution(cell->getCellCenter(), 0) << endl;

//    output << "SCALARS e double" << endl;
//    output << "LOOKUP_TABLE default" << endl;

//    for (const shared_ptr<Cell> cell : mesh.cells)
//        output << cell->reconstructSolution(cell->getCellCenter(), 4) << endl;

//    output << "SCALARS p double" << endl;
//    output << "LOOKUP_TABLE default" << endl;

//    for (const shared_ptr<Cell> cell : mesh.cells)
//        output << problem.getPressure(cell->reconstructSolution(cell->getCellCenter())) << endl;


//    output << "VECTORS U double" << endl;

//    for (const shared_ptr<Cell> cell : mesh.cells)
//    {
//        double rho = cell->reconstructSolution(cell->getCellCenter(), 0);
//        output << cell->reconstructSolution(cell->getCellCenter(), 1) / rho << endl;
//        output << cell->reconstructSolution(cell->getCellCenter(), 2) / rho << endl;
//        output << 0.0 << endl;
//    }


    // get point data

    output << "POINT_DATA " << mesh.nEntitiesTotal << endl;

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


    output.close();
}




void Solver::setInitialConditions()
{
    int nCells = mesh.nCells;

    problem.alpha.resize(nCells);

    numvector<double, 5*nShapes> rhs;

//#pragma omp parallel for
    for (int k = 0; k < nCells; ++k)
    {
        rhs = mesh.cells[k]->projection(problem.init);
        problem.alpha[k] = mesh.cells[k]->correctNonOrtho(rhs);
        alphaPrev[k] = problem.alpha[k];
    }

} // end setInitialConditions

void Solver::setDefinedCoefficients(string fileName)
{
    ifstream reader;
    reader.open(fileName);

    if (!reader.is_open())
    {
        cout << "File " << fileName << " is not found\n";
        exit(0);
    }

    int nCells = mesh.nCells;

    problem.alpha.resize(nCells);

    numvector<double, 5*nShapes> rhs;

    for (int k = 0; k < nCells; ++k)
    {
        for (int j = 0; j < 5*nShapes; ++j)
            reader >> problem.alpha[k][j];

        alphaPrev[k] = problem.alpha[k];
    }

    reader.close();
}

void Solver::setMeshPointerForDiagBC()
{
//    int nEdgesBou = mesh.edgesBoundary.size();

//    for (int i = 0; i < nEdgesBou; ++i)
//    {
//        shared_ptr<EdgeBoundaryDiagProjection> edgeBound = dynamic_pointer_cast<EdgeBoundaryDiagProjection>(mesh.edgesBoundary[i]);

//        if (edgeBound == NULL)
//            continue;

//        edgeBound->setMeshPointer(mesh);
//    }
} // end setMeshPointerForDiagBC

vector<numvector<double, 5 * nShapes>> Solver::assembleRHS(const std::vector<numvector<double, 5 * nShapes> > &alpha)
{
    problem.setAlpha(alpha);

    //int nEdgesInt = mesh.edgesInternal.size();
    //int nEdgesBou = mesh.edgesBoundary.size();

    int nCells = mesh.cells.size();

    // compute fluxes in gauss points on edges
#pragma omp parallel for 
    for (size_t i = 0; i < mesh.edges.size(); ++i)
        mesh.edges[i]->getLocalFluxes(flux);


    vector<numvector<double, 5 * nShapes>> rhs(mesh.nCells);

#pragma omp parallel for
    for (int k = 0; k < nCells; ++k) // for all cells
    {
        // compute internal integral
        rhs[k] = mesh.cells[k]->cellIntegral();

        // compute boundary integrals
        for (int i = 0; i < mesh.cells[k]->nEntities; ++i)
            rhs[k] -= mesh.cells[k]->edges[i]->boundaryIntegral(mesh.cells[k]);
    }

    return rhs;
}

vector<numvector<double, 5 * nShapes>> Solver::correctNonOrtho(std::vector<numvector<double, 5 * nShapes>> &alpha) const
{
    vector<numvector<double, 5 * nShapes>> res = alpha;

#pragma omp parallel for
    for (int i = 0; i < mesh.nCells; ++i)
        res[i] = mesh.cells[i]->correctNonOrtho(alpha[i]);

    return res;
}

vector<numvector<double, 5 * nShapes>> Solver::correctPrevIter(std::vector<numvector<double, 5 * nShapes>> &alpha) const
{
    vector<numvector<double, 5 * nShapes>> res = alpha;

#pragma omp parallel for
    for (int i = 0; i < mesh.nCells; ++i)
        res[i] = mesh.cells[i]->correctPrevIter(alpha[i]);

    return res;
}
