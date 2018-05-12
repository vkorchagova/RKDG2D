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
    output.open(fileName + ".vtk");

    mesh.exportMeshVTK(output);

    output << "CELL_DATA " << mesh.nCells << endl;

    output << "SCALARS rho double" << endl;
    output << "LOOKUP_TABLE default" << endl;
    //output << "rho" << " 1" << " " << mesh.nCells << " float" << endl;

    for (int i = 0; i < mesh.nCells; ++i)
       output << mesh.cells[i]->reconstructSolution(mesh.cells[i]->getCellCenter(),0) << endl;

    output << "SCALARS e double" << endl;
    output << "LOOKUP_TABLE default" << endl;
    //output << "e" << " 1" << " " << mesh.nCells << " float" << endl;

    for (int i = 0; i < mesh.nCells; ++i)
       output << mesh.cells[i]->reconstructSolution(mesh.cells[i]->getCellCenter(),4) << endl;

    output << "VECTORS rhoU double" << endl;
    //output << "rhoU " << " 3" << " " << mesh.nCells << " float" << endl;

    for (int i = 0; i < mesh.nCells; ++i)
    {
       output << mesh.cells[i]->reconstructSolution(mesh.cells[i]->getCellCenter(),1) << ' ';
       output << mesh.cells[i]->reconstructSolution(mesh.cells[i]->getCellCenter(),2) << ' ';
       output << 0.0 << endl;
    }

    output << "POINT_DATA " << mesh.nodes.size() << endl;

    output.close();

}


void Solver::setInitialConditions()
{
    int nCells = mesh.nCells;

    problem.alpha.resize(nCells);

    numvector<double, 5*nShapes> rhs;

    for (int k = 0; k < nCells; ++k)
    {
        rhs = mesh.cells[k]->projection(problem.init);
        problem.alpha[k] = mesh.cells[k]->correctNonOrtho(rhs);
        alphaPrev[k] = problem.alpha[k];
    }

//    ofstream writer;
//    writer.open("alphaCoeffs/0.000000");

//    write(writer, problem.alpha);

//    writer.close();

} // end setInitialConditions

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

    for (size_t i = 0; i < mesh.edges.size(); ++i)
        mesh.edges[i]->getLocalFluxes(flux);


    vector<numvector<double, 5 * nShapes>> rhs(mesh.nCells);

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

    for (int i = 0; i < mesh.nCells; ++i)
        res[i] = mesh.cells[i]->correctNonOrtho(alpha[i]);

    return res;
}

vector<numvector<double, 5 * nShapes>> Solver::correctPrevIter(std::vector<numvector<double, 5 * nShapes>> &alpha) const
{
    vector<numvector<double, 5 * nShapes>> res = alpha;

    for (int i = 0; i < mesh.nCells; ++i)
        res[i] = mesh.cells[i]->correctPrevIter(alpha[i]);

    return res;
}
