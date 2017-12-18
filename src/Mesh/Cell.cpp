#include "Cell.h"
#include "Edge.h"
#include <iostream>

using namespace std;

// ------------------ Constructors & Destructor ----------------

Cell::Cell(const numvector<std::shared_ptr<Edge>, nEdges> &defEdges)
{
    edges = defEdges;
    
    step.x() = edges[0]->nodes[1]->x() - edges[0]->nodes[0]->x();
    step.y() = edges[nEdges-1]->nodes[1]->y() - edges[nEdges-1]->nodes[0]->y();
    
    area = step.x() * step.y();

    J = 0.25 * area;
    
    center = *(edges[0]->nodes[0]) + 0.5 * step;
    
    setBasisFunctions();
    setGaussPoints();
}

// ------------------ Private class methods --------------------


bool Cell::insideCell(const Point& point) const
{
    double epsilon = 1e-10;

    bool xCond = (center.x() - 0.5*step.x() - epsilon) < point.x()  && (center.x() + 0.5*step.x() + epsilon) > point.x();
    bool yCond = (center.y() - 0.5*step.y() - epsilon) < point.y()  && (center.y() + 0.5*step.y() + epsilon) > point.y();

    if (xCond && yCond)
        return true;

    return false;
}

Point Cell::localToGlobal(const Point& point) const
{
    return Point( { 0.5 * step.x() * point.x() + center.x(), 0.5 * step.y() * point.y() + center.y() } );
}

void Cell::setGaussPoints()
{
    double isqrt3 = 1.0/1.732050807568877;

    gPoints2D[0] = localToGlobal(Point({-isqrt3, -isqrt3}));
    gPoints2D[1] = localToGlobal(Point({ isqrt3, -isqrt3}));
    gPoints2D[2] = localToGlobal(Point({-isqrt3,  isqrt3}));
    gPoints2D[3] = localToGlobal(Point({ isqrt3,  isqrt3}));

    gWeights2D = { 1.0, 1.0, 1.0, 1.0 };
}

void Cell::setBasisFunctions()
{
    double& hx = step.x();
    double& hy = step.y();
    
    phi.emplace_back([&](const Point& r){ return 1.0 / sqrt(hx*hy); });
    phi.emplace_back([&](const Point& r){ return sqrt(12.0 / hy / pow(hx, 3)) * (r.x() - center.x()); });
    phi.emplace_back([&](const Point& r){ return sqrt(12.0 / hx / pow(hy, 3)) * (r.y() - center.y()); });

    gradPhi.emplace_back([&](const Point& r)->Point { return Point({0.0, 0.0}); });
    gradPhi.emplace_back([&](const Point& r)->Point { return Point({sqrt(12.0 / hy / pow(hx, 3)), 0.0}); });
    gradPhi.emplace_back([&](const Point& r)->Point { return Point({0.0, sqrt(12.0 / hx / pow(hy, 3))}); });
}

// ------------------ Public class methods ---------------------

// ----- geometric -----

vector<shared_ptr<Point>> Cell::getCellCoordinates() const
{
    vector<shared_ptr<Point>> nodeCoordinates(nEdges);

    for (int i = 0; i < nEdges; ++i)
    {
        nodeCoordinates[i] = edges[i]->nodes[0];
    }

    return nodeCoordinates;

} // end getCellCoordinates


// ----- RKDG -----

void Cell::setProblem(const Problem& prb)
{
    problem = &prb;

} // end setProblem

numvector<double, 5> Cell::reconstructSolution(const Point& point) const
{
    if (!insideCell(point))
    {
        std::cout << "Error: point (" << point.x() << ", " << point.y() << ") is not inside cell #" << number << std::endl;
        std::cout << "Cell nodes:" << std::endl;

        vector<shared_ptr<Point>> ccoord = getCellCoordinates();

        for (int i = 0; i < ccoord.size(); ++i)
            std::cout << "(" << ccoord[i]->x() << "; " << ccoord[i]->y() << ")" << endl;

        exit(1);
    }

    numvector<double, 5> sol(0.0);

    for (int i = 0; i < 5; ++i)
    {
        for (int j = 0; j < nShapes; ++j)
            sol[i] += phi[j](point) * problem->alpha[number][i * nShapes + j];
    }

    return sol;

} // end reconstructSolution


numvector<double, 5 * nShapes> Cell::getLocalInitialConditions(std::function<numvector<double,5>(const Point& point)>& init) const
{
    numvector<double, 5 * nShapes> alpha;

    for (int q = 0; q < nShapes; ++q)
    {
        std::function<numvector<double, 5>(const Point&)> f = \
                [&](const Point& p) {  return phi[q](p) * init(p); };

        numvector<double, 5> buffer = integrate(f);

        for (int p = 0; p < 5; ++p)
        {
            alpha[p*nShapes + q] = buffer[p];
        }// for p
    }

    return alpha;

} // end setLocalInitialConditions

numvector<double, 5> Cell::integrate( const std::function<numvector<double, 5>(const Point &)> &f) const
{
    numvector<double, 5> res = 0.0;

    double J = area*0.25;

    for (int i = 0; i < nGP; ++i)
    {
        numvector<double,5> resF = f(gPoints2D[i]);

        for (int k = 0; k < 5; ++k)
        {
           res[k] += gWeights2D[i] * resF[k];
        }
    }

    return res * J;

} // end integrate 2D of vector function

numvector<double, 5 * nShapes> Cell::cellIntegral()
{
    numvector<double, 5> sol;
    numvector<double, 5> resV;
    numvector<double, 5 * nShapes> res (0.0);

    for (int i = 0; i < nGP; ++i)
    {
        sol = reconstructSolution(gPoints2D[i]);

        for (int q = 0; q < nShapes; ++q)
        {
            resV = problem->fluxF(sol) * gradPhi[q](gPoints2D[i])[0] + \
                   problem->fluxG(sol) * gradPhi[q](gPoints2D[i])[1];

            for (int p = 0; p < 5; ++p)
                res[p * nShapes + q] += resV[p] * gWeights2D[i];
        }
    }

    return res * J;
}

