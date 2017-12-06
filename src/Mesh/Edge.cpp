#include "Edge.h"

using namespace std;

Edge::Edge()
{

}

Edge::~Edge()
{

}

Edge::Edge(const Edge& rhs)
{
    nodes = rhs.nodes;
    neibCells = rhs.neibCells;
    localFluxes = rhs.localFluxes;
}

Edge& Edge::operator=(const Edge& rhs)
{
    nodes = rhs.nodes;
    neibCells = rhs.neibCells;
    localFluxes = rhs.localFluxes;

    return *this;
}
