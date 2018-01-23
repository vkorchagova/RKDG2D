#include "Solver.h"
#include <iostream>

using namespace std;

Solver::Solver( Mesh2D& msh, Problem &prb, Flux& flx) : mesh(msh),problem(prb),flux(flx)
{
    alphaPrev.resize(mesh.nCells);
    alphaNext.resize(mesh.nCells);
}


void Solver::write(ostream& writer, const vector<numvector<double,5*nShapes>>& coeffs) const
{
    writer.precision(12);

    for (size_t k = 0; k < coeffs.size(); ++k)
    {
        for (int i = 0; i < 5*nShapes; ++i)
            writer << coeffs[k][i] << ' ';

        writer << endl;
    }
} //end write

void Solver::setInitialConditions()
{
    int nCells = mesh.nCells;

    problem.alpha.resize(nCells);

    for (int k = 0; k < nCells; ++k)
    {
        mesh.cells[k]->setProblem(problem);
        problem.alpha[k] = mesh.cells[k]->projection(problem.init);
        alphaPrev[k] = problem.alpha[k];
    }

    ofstream writer;
    writer.open("alphaCoeffs/0.000000");

    write(writer, problem.alpha);

    writer.close();

} // end setInitialConditions

void Solver::setBoundaryConditions() const
{
    int nEdgesHor = mesh.edgesHor.size();
    int nEdgesVer = mesh.edgesVer.size();

    for (int i = 0; i < nEdgesHor; ++i)
        mesh.edgesHor[i]->setBoundaryFunction(problem.infty);

    for (int i = 0; i < nEdgesVer; ++i)
        mesh.edgesVer[i]->setBoundaryFunction(problem.infty);
} // setBoundaryConditions

void Solver::setMeshPointerForDiagBC()
{
    int nEdgesHor = mesh.edgesHor.size();
    int nEdgesVer = mesh.edgesVer.size();

    for (int i = 0; i < nEdgesHor; ++i)
    {
        shared_ptr<EdgeBoundaryDiagProjection> edgeBound = dynamic_pointer_cast<EdgeBoundaryDiagProjection>(mesh.edgesHor[i]);

        if (edgeBound == NULL)
            continue;

        edgeBound->setMeshPointer(mesh);
    }

    for (int i = 0; i < nEdgesVer; ++i)
    {
        shared_ptr<EdgeBoundaryDiagProjection> edgeBound = dynamic_pointer_cast<EdgeBoundaryDiagProjection>(mesh.edgesVer[i]);

        if (edgeBound == NULL)
            continue;

        edgeBound->setMeshPointer(mesh);
    }
} // end setMeshPointerForDiagBC

vector<numvector<double, 5 * nShapes>> Solver::assembleRHS(const std::vector<numvector<double, 5 * nShapes> > &alpha)
{
    problem.setAlpha(alpha);

    int nEdgesHor = mesh.edgesHor.size();
    int nEdgesVer = mesh.edgesVer.size();

    int nCells = mesh.cells.size();

    // compute fluxes in gauss points on edges

    for (int i = 0; i < nEdgesHor; ++i)
        mesh.edgesHor[i]->getLocalFluxes(flux);

    for (int i = 0; i < nEdgesVer; ++i)
        mesh.edgesVer[i]->getLocalFluxes(flux);


    vector<numvector<double, 5 * nShapes>> rhs(mesh.nCells);

    for (int k = 0; k < nCells; ++k) // for all cells
    {
        // compute internal integral
        rhs[k] = mesh.cells[k]->cellIntegral();

        // compute boundary integrals
        for (int i = 0; i < 4; ++i)
            rhs[k] -= mesh.cells[k]->edges[i]->boundaryIntegral(mesh.cells[k]);
    }

    return rhs;
}
