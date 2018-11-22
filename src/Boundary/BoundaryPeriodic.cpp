#include "BoundaryPeriodic.h"
#include "Edge.h"

numvector<double, 5> BoundaryPeriodic::applyBoundary(const numvector<double, 5>& solLeft, const Point& n, int numGP) const
{
    return adjEdge->neibCells[0]->reconstructSolution(adjEdge->gPoints[numGP]);
}
