#include "Edge.h"
#include <iostream>

using namespace std;

// ------------------ Constructors & Destructor ----------------

Edge::Edge(const Point& p1, const Point& p2)
{
    nodes[0] = make_shared<Point>(p1);
    nodes[1] = make_shared<Point>(p2);

    const double isqrt3 = 1.0/1.732050807568877;

    Point c = 0.5 * (p2 + p1);
    Point m = 0.5 * (p2 - p1);

    // gauss points in global coordinate system
    gPoints[0] = c - isqrt3 * m;
    gPoints[1] = c + isqrt3 * m;
    
    // weights for gauss integration
    gWeights = { 1.0, 1.0 };
    
    length = (p2 - p1).length();

    //- jacobian
    J = 0.5 * length ;
}

// ------------------ Private class methods --------------------


// ------------------ Public class methods ---------------------

//// RKDG methods
numvector<double, 5 * nShapes> Edge::boundaryIntegral(const std::shared_ptr<Cell> &cell) const
{
    numvector<double, 5 * nShapes> res (0.0);

    double sign = (cell == neibCells[0]) ? 1.0 : -1.0;

    for (int i = 0; i < nGP; ++i)
        for (int q = 0; q < nShapes; ++q)
            for (int p = 0; p < 5; ++p)
                res[p*nShapes + q] += localFluxes[i][p] * ( gWeights[i] * cell->phi[q](gPoints[i]) );

    return res * J * sign;
}

