#include "Cell.h"
#include <iostream>

// ------------------ Constructors & Destructor ----------------

Cell::Cell()
{

}

Cell::Cell(numvector<Edge*, 4>& defEdges)
{
    edges = defEdges;

    getSteps();
    getArea();
    getCellCenter();

}

Cell::~Cell()
{

}

// ------------------ Private class methods --------------------

void Cell::getSteps()
{
    hx = edges[0]->nodes[1]->x() - edges[0]->nodes[0]->x();
    hy = edges[nEdges-1]->nodes[1]->y() - edges[nEdges-1]->nodes[0]->y();
}

void Cell::getArea()
{
    area = hx * hy;
}

void Cell::getCellCenter()
{
    Point c (edges[0]->nodes[0]->x() + 0.5*hx, edges[0]->nodes[0]->y() + 0.5*hy);

    center = c;
}

// ------------------ Public class methods ---------------------

numvector<Point*, 4> Cell::getCellCoordinates()
{
    numvector<Point*, 4> nodeCoordinates;

    for (int i = 0; i < nEdges; ++i)
    {
        nodeCoordinates[i] = edges[i]->nodes[0];
    }

    return nodeCoordinates;
}
