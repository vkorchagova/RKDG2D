#include "Edge.h"
#include <iostream>

using namespace std;

// ------------------ Constructors & Destructor ----------------

Edge::Edge(const Node &p1, const Node &p2)
{
    nodes[0] = make_shared<Node>(p1);
    nodes[1] = make_shared<Node>(p2);

    Point c = 0.5 * (p2 + p1);
    Point m = 0.5 * (p2 - p1);

    // gauss points in global coordinate system & weights
    if (nShapes == 6)
    {
        nGP = 3;

        const double sqrtfrac35 = 0.7745966692414834;

        gPoints.push_back( c - sqrtfrac35 * m );
        gPoints.push_back( c );
        gPoints.push_back( c + sqrtfrac35 * m );

        const double frac59 = 0.5555555555555556;
        const double frac89 = 0.8888888888888889;

        gWeights = { frac59, frac89, frac59 };

    }
    else if (nShapes == 3 || nShapes == 1)
    {
        // pure Gauss
        nGP = 2;

        const double isqrt3 = 0.57735026918962576;

        gPoints.push_back( c - isqrt3 * m );
        gPoints.push_back( c + isqrt3 * m );

        // weights for gauss integration
        gWeights = { 1.0, 1.0 };

        // Gauss --- Lobatto
//        nGP = 3;

//        gPoints.push_back( c - m );
//        gPoints.push_back( c );
//        gPoints.push_back( c + m );

//        gWeights = { 1.0/3.0, 4.0/3.0, 1.0/3.0 };

    }
    else
    {
        cout << "Wrong number of basis functions " << nShapes;
        cout << "Avaliable: 1,3,6 FF in 2D case \n";

        exit(0);
    }

    localFluxes.resize(nGP);
    
    length = (p2 - p1).length();

    /// jacobian
    J = 0.5 * length ;
}

// ------------------ Private class methods --------------------




// ------------------ Public class methods ---------------------

//// RKDG methods

numvector<double, 5 * nShapes> Edge::boundaryIntegral(const std::shared_ptr<Cell> &cell) const
{
    numvector<double, 5 * nShapes> res (0.0);
    double gW = 0.0;

    double sign = (cell == neibCells[0]) ? 1.0 : -1.0;

    for (int i = 0; i < nGP; ++i)
    {
        gW = gWeights[i];
        
        for (int q = 0; q < nShapes; ++q)
            for (int p = 0; p < 5; ++p)
                res[p*nShapes + q] += localFluxes[i][p] * ( gW * cell->phi[q](gPoints[i]) );
    
    }

    return res * J * sign;
}


double Edge::getMassFlux(const std::shared_ptr<Cell> &cell) const
{
    double res = 0.0;

    double sign = (cell == neibCells[0]) ? 1.0 : -1.0;

    for (int i = 0; i < nGP; ++i)
        res += localFluxes[i][0] * gWeights[i] / cell->reconstructSolution(gPoints[i],0);

    return res * J * sign;
}


