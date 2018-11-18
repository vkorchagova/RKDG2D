#include "EdgeBoundaryDiagProjection.h"

using namespace std;

// ------------------ Constructors & Destructor ----------------

// ------------------ Private class methods ---------------------

Point EdgeBoundaryDiagProjection::getDiagProjection(const Point& p) const
{
    Point diag = Point({0,0});// ({ neibCells[0]->h().x(), neibCells[0]->h().y()});
    diag.normalize();

    double projection = diag * p;

    return diag * projection;

}

int EdgeBoundaryDiagProjection::getNumDiagCell(const Point& p)
{
    for (shared_ptr<Cell> cell : mesh->cells)
    {
        if (cell->insideCell(p))
            return cell->number;
    }

    cout << "No cell" << endl;
    exit(0);
}

// ------------------ Public class methods ---------------------

void EdgeBoundaryDiagProjection::getLocalFluxes(const Flux &flux)
{

    for (int i = 0; i < nGP; ++i)
    {
        numvector<double, 5> solInner = rotate(neibCells[0]->reconstructSolution(gPoints[i]),n);

        Point proj = getDiagProjection(gPoints[i]);
        int numCell = getNumDiagCell(proj);
        numvector<double, 5> solOuter = rotate(mesh->cells[numCell]->reconstructSolution(proj),n);

        localFluxes[i] = flux.evaluate(solInner, solOuter, n);
    }
}
