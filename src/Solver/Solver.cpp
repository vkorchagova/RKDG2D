#include "Solver.h"

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

    for (int k = 0; k < nCells; ++k)
    {
        mesh->cells[k].setProblem(*problem);
        mesh->cells[k].setLocalInitialConditions(problem->init);

        write(writer,mesh->cells[k].alphaPrev);
    }

    writer.close();

} // end setInitialConditions

void Solver::assembleRHS()
{
    int nEdgesHor = mesh->edgesHor.size();
    int nEdgesVer = mesh->edgesVer.size();

    int nCells = mesh->cells.size();

    for (int i = 0; i < nEdgesHor; ++i)
        mesh->edgesHor[i].getLocalFluxes();

    numvector<double, 5*nShapes> res;


    for (int i = 0; i < nCells; ++i) // for all cells
    {
        res = mesh->cells[i].getLocalRHS();
    }





}
