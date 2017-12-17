#include "Solver.h"
#include <iostream>

using namespace std;

Solver::Solver( Mesh2D& msh, Problem &prb, Flux& flx) : mesh(msh),problem(prb),flux(flx)
{
    alphaPrev.resize(mesh.nCells);
    alphaNext.resize(mesh.nCells);
}


void Solver::write(ostream& writer, const numvector<double,5*nShapes>& coeffs) const
{
    for (int i = 0; i < 5*nShapes; ++i)
        writer << coeffs[i] << ' ';

    writer << endl;
} //end write

void Solver::setInitialConditions()
{
    int nCells = mesh.nCells;

    ofstream writer;
    writer.open("alphaCoeffs/0.000000");

    problem.alpha.resize(nCells);

    for (int k = 0; k < nCells; ++k)
    {
        mesh.cells[k]->setProblem(problem);
        problem.alpha[k] = mesh.cells[k]->getLocalInitialConditions(problem.init);
        alphaPrev[k] = problem.alpha[k];

        write(writer, problem.alpha[k]);
    }

    writer.close();

} // end setInitialConditions

void Solver::initBoundaryConditions() const
{
    int nEdgesHor = mesh.edgesHor.size();
    int nEdgesVer = mesh.edgesVer.size();

    for (int i = 0; i < nEdgesHor; ++i)
        mesh.edgesHor[i]->setBoundaryFunction(problem.infty);

    for (int i = 0; i < nEdgesVer; ++i)
        mesh.edgesVer[i]->setBoundaryFunction(problem.infty);
}

void Solver::assembleRHS(const std::vector<numvector<double, 5 * nShapes> > &alpha)
{
    problem.getAlpha(alpha);

    int nEdgesHor = mesh.edgesHor.size();
    int nEdgesVer = mesh.edgesVer.size();

    int nCells = mesh.cells.size();

    // compute fluxes in gauss points on edges

    for (int i = 0; i < nEdgesHor; ++i)
        mesh.edgesHor[i]->getLocalFluxesHor(flux);

    for (int i = 0; i < nEdgesVer; ++i)
        mesh.edgesVer[i]->getLocalFluxesVer(flux);

    for (int i = 0; i < nEdgesHor; ++i)
        cout << "fluxes hor: " <<  mesh.edgesHor[i]->localFluxes << endl;

    for (int i = 0; i < nEdgesVer; ++i)
        cout << "fluxes ver: " <<  mesh.edgesVer[i]->localFluxes << endl;

    // compute boundary integrals

    numvector<double, 5 * nShapes> res;

    for (int k = 0; k < nCells; ++k) // for all cells
    {
        cout << "Cell #" << k << endl;

        alphaNext =

        for (int i = 0; i < 4; ++i)
            alphaNext[k] += mesh.cells[k]->edges[i]->boundaryIntegral(mesh.cells[k]);



        cout << alphaNext[k] << endl;
    }


}
