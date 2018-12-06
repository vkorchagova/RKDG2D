#include "Edge.h"
#include <iostream>

using namespace std;

// ------------------ Constructors & Destructor ----------------

Edge::Edge(const Point &p1, const Point &p2)
//Edge::Edge(const vector<shared_ptr<Point>> &p) : nodes(p)
{
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
   
	// ...
    length = (p2 - p1).length();

    //- jacobian
    J = 0.5 * length;
}

