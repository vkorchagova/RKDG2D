#include "Edge.h"
#include <iostream>

using namespace std;

// ------------------ Constructors & Destructor ----------------

Edge::Edge(const Point& p1, const Point& p2, const Flux& flux_) : flux(flux_)
{
    nodes[0] = &p1;
    nodes[1] = &p2;

    const double isqrt3 = 1.0/1.732050807568877;

    Point c = 0.5 * (p2 + p1);
    Point m = 0.5 * (p2 - p1);

    // gauss points in global coordinate system
    gPoints[0] = c - isqrt3 * m;
    gPoints[1] = c + isqrt3 * m;
    
    // weights for gauss integration
    gWeights = { 1.0, 1.0 };
}

// ------------------ Private class methods --------------------


// ------------------ Public class methods ---------------------

//// RKDG methods
numvector<double, 5> Edge::boundaryIntegral(const std::function<double(const Point&)>& phi) const
{
    numvector<double, 5> res (0.0);

//    std::cout << "boundary integral: ";

    for (int i = 0; i < nGP; ++i)
    {
        for (int k = 0; k < 5; ++k)
        {
            res[k] += ( gWeights[i] * phi(gPoints[i]) ) * localFluxes[i][k];
        }
    }

    return res;
}

