#include "Solver.h"
#include <iostream>

using namespace std;

void Solver::write(ostream& writer, const numvector<double,5*nShapes>& coeffs)
{
    for (int i = 0; i < 5*nShapes; ++i)
        writer << coeffs[i] << ' ';

    writer << endl;
} //end write

void Solver::setInitialConditions()
{
    int nCells = mesh->nCells;

    ofstream writer;
    writer.open("alphaCoeffs/0.000000");

    problem->alpha.resize(nCells);

    for (int k = 0; k < nCells; ++k)
    {
        mesh->cells[k].setProblem(*problem);
        problem->alpha[k] = mesh->cells[k].setLocalInitialConditions(problem->init);

        write(writer,problem->alpha[k]);
    }

    writer.close();

} // end setInitialConditions

void Solver::initBoundaryConditions()
{
    int nBE = mesh->edgesBound.size();

    for (int i = 0; i < nBE; ++i)
        mesh->edgesBound[i].infty = problem->infty;
}

void Solver::initFluxes(Flux& flux)
{
    int nEdgesHor = mesh->edgesHor.size();
    int nEdgesVer = mesh->edgesVer.size();

    for (int i = 0; i < nEdgesHor; ++i)
        mesh->edgesHor[i].setFlux(flux);

    for (int i = 0; i < nEdgesVer; ++i)
        mesh->edgesVer[i].setFlux(flux);
}

void Solver::assembleRHS(std::vector<numvector<double, 5 * nShapes> > &alpha)
{
    problem->getAlpha(alpha);

    int nEdgesHor = mesh->edgesHor.size();
    int nEdgesVer = mesh->edgesVer.size();

    int nCells = mesh->cells.size();

    // compute fluxes in gauss points on edges

    for (int i = 0; i < nEdgesHor; ++i)
    {
        mesh->edgesHor[i].getLocalFluxes();
        //std::cout << mesh->edgesHor[i].localFluxes[0] << '\n';
    }


    for (int i = 0; i < nEdgesVer; ++i)
        mesh->edgesVer[i].getLocalFluxes();

    // compute boundary integrals

    numvector<double, 5 * nShapes> res;

    for (int i = 0; i < nCells; ++i) // for all cells
    {
        res = mesh->cells[i].getLocalRHS();
    }


}
