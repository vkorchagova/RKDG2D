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
    setBasisFunctions();
    setGaussPoints();
}

Cell::~Cell()
{

}

Cell::Cell(const Cell& c)
{
    edges = c.edges;
    area = c.area;
    center = c.center;
    hx = c.hx;
    hy = c.hy;
    number = c.number;

    setBasisFunctions();
    setGaussPoints();
}

Cell& Cell::operator=(const Cell& c)
{
    edges = c.edges;
    area = c.area;
    center = c.center;
    hx = c.hx;
    hy = c.hy;
    number = c.number;

    setBasisFunctions();
    setGaussPoints();

    return *this;
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

bool Cell::insideCell(Point& point)
{
    bool xCond = (center.x() - 0.5*hx) < point.x() && (center.x() + 0.5*hx) > point.x();
    bool yCond = (center.y() - 0.5*hy) < point.y() && (center.y() + 0.5*hy) > point.y();

    if (xCond && yCond)
        return true;

    return false;
}

Point Cell::localToGlobal(Point &point)
{
    Point global (0.5 * hx * point.x() + center.x(), 0.5*hy*point.y() + center.y());

    return global;
}

void Cell::setGaussPoints()
{

    nGP = 4;

    double sqrt3 = 1.0/1.732050807568877;

    gPoints2D[0].set(-sqrt3, -sqrt3);
    gPoints2D[1].set( sqrt3, -sqrt3);
    gPoints2D[2].set(-sqrt3,  sqrt3);
    gPoints2D[3].set( sqrt3,  sqrt3);

    gWeights2D = { 1.0, 1.0, 1.0, 1.0 };
}

void Cell::setBasisFunctions()
{
    phi.resize(3);

    phi[0] = [=](const Point& r){ return 1.0 / sqrt(hx*hy); };
    phi[1] = [=](const Point& r){ return sqrt(12.0 / hy / pow(hx, 3)) * (r.x() - center.x()); };
    phi[2] = [=](const Point& r){ return sqrt(12.0 / hx / pow(hy, 3)) * (r.y() - center.y()); };

    gradPhi[0] = [=](const Point& r){ Point p (0.0, 0.0);                          return p; };
    gradPhi[1] = [=](const Point& r){ Point p (sqrt(12.0 / hy / pow(hx, 3)), 0.0); return p; };
    gradPhi[2] = [=](const Point& r){ Point p (0.0, sqrt(12.0 / hx / pow(hy, 3))); return p; };
}

// ------------------ Public class methods ---------------------

// ----- geometric -----

numvector<Point*, 4> Cell::getCellCoordinates()
{
    numvector<Point*, 4> nodeCoordinates;

    for (int i = 0; i < nEdges; ++i)
    {
        nodeCoordinates[i] = edges[i]->nodes[0];
    }

    return nodeCoordinates;

} // end getCellCoordinates


// ----- RKDG -----

void Cell::setProblem(Problem& prb)
{
    problem = &prb;

} // end setProblem

numvector<double, 5> Cell::reconstructSolution(const numvector<double, nShapes * 5>& alpha, Point& point)
{
    if (!insideCell(point))
    {
        std::cout << "Error: point (" << point.x() << ", " << point.y() << "not inside cell #" << number << std::endl;
        exit(1);
    }

    numvector<double, 5> sol(0.0);

    for (int i = 0; i < 5; ++i)
    {
        for (int j = 0; j < nShapes; ++j)
            sol[i] += phi[j](point)*alpha[i*nShapes + j];
    }

    return sol;

} // end reconstructSolution


void Cell::setLocalInitialConditions(std::function<numvector<double,5>(const Point& point)>& init)
{
    for (int q = 0; q < nShapes; ++q)
    {
        std::function<numvector<double, 5>(const Point&)> f = \
                [=](const Point& p) {  return phi[q](p) * init(p); };

        numvector<double, 5> buffer = integrate(f);

        for (int p = 0; p < 5; ++p)
        {
            alphaPrev[p*nShapes + q] = buffer[p];
        }// for p
    }

} // end setLocalInitialConditions

numvector<double,5> Cell::integrate( const std::function<numvector<double, 5>(Point &)> &f)
{
    numvector<double,5> res = 0.0;

    double J = area*0.25;

    for (int i = 0; i < nGP; ++i)
    {
        Point glob = localToGlobal(gPoints2D[i]);

        numvector<double,5> resF = f(glob);

        for (int k = 0; k < 5; ++k)
        {
           res[k] += gWeights2D[i] * resF[k];
        }
    }

    return res*J;

} // end integrate 2D of vector function
