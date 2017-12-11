#include "Edge.h"
#include <iostream>

using namespace std;

// ------------------ Constructors & Destructor ----------------

Edge::Edge(Point* p1, Point* p2)
{
    nodes[0] = p1;
    nodes[1] = p2;

    double sqrt3 = 1.0/1.732050807568877;

    double xc = 0.5 * (p2->x() + p1->x());
    double yc = 0.5 * (p2->y() + p1->y());
    double xm = 0.5 * (p2->x() - p1->x());
    double ym = 0.5 * (p2->y() - p1->y());

    // gauss points in global coordinate system
    gPoints[0].set( - sqrt3 * xm + xc, - sqrt3 * ym + yc );
    gPoints[1].set(   sqrt3 * xm + xc,   sqrt3 * ym + yc );

    // weights for gauss integration
    gWeights = { 1.0, 1.0 };
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

// ------------------ Private class methods --------------------



// ------------------ Public class methods ---------------------

//// RKDG methods

void Edge::setFlux(Flux& flx)
{
    flux = &flx;
}

numvector<double, 5> Edge::boundaryIntegral(std::function<double(const Point&)>& phi)
{
    numvector<double, 5> res (0.0);

//    std::cout << "boundary integral: ";

    for (int i = 0; i < nGP; ++i)
    {
        for (int k = 0; k < 5; ++k)
        {
            res[k] += ( gWeights[i] * phi(gPoints[i]) ) * localFluxes[i][k] ;
        }
    }

//    for (int k = 0; k < 5; ++k)
//        std::cout << res[k] << ' ';
//    std::cout << std::endl;

    return res;
}

