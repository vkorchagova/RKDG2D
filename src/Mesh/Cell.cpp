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
    double epsilon = 1e-16;

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
    if (nShapes == 6)
    {
        nGP = 9;

        const double sqrtfrac35 = 0.7745966692414834;

        for (int i = -1; i <= 1; i++)
            for (int j = -1; j <= 1; j++)
                gPoints2D.push_back(localToGlobal(Point({ i * sqrtfrac35, j * sqrtfrac35 })));

        const double frac59 = 0.5555555555555556;
        const double frac89 = 0.8888888888888889;

        vector <double> gWeights1D = { frac59, frac89, frac59 };

        for (int i = 0; i <= 2; ++i)
            for (int j = 0; j <= 2; ++j)
                gWeights2D.push_back(gWeights1D[i] * gWeights1D[j]);

    }
    else if (nShapes == 3 || nShapes == 1)
    {
        nGP = 4;

        double isqrt3 = 0.57735026918962576;

        for (int i = -1; i <= 1; i += 2)
            for (int j = -1; j <= 1; j += 2)
                gPoints2D.push_back(localToGlobal(Point({ i * isqrt3, j * isqrt3 })));

        gWeights2D = { 1.0, 1.0, 1.0, 1.0 };
    }
    else
    {
        cout << "Wrong number of basis functions " << nShapes;
        cout << "Avaliable: 1,3,6 FF in 2D case \n";

        exit(0);
    }
}

void Cell::setBasisFunctions()
{
    double& hx = step.x();
    double& hy = step.y();

    offsetPhi.push_back(1.0 / sqrt(hx*hy));
    offsetPhi.push_back(sqrt(12.0 / hy / pow(hx, 3)));
    offsetPhi.push_back(sqrt(12.0 / hx / pow(hy, 3)));
    
    phi.emplace_back([&](const Point& r){ return 1.0 / sqrt(hx*hy); });
    phi.emplace_back([&](const Point& r){ return sqrt(12.0 / hy / pow(hx, 3)) * (r.x() - center.x()); });
    phi.emplace_back([&](const Point& r){ return sqrt(12.0 / hx / pow(hy, 3)) * (r.y() - center.y()); });

    gradPhi.emplace_back([&](const Point& r)->Point { return Point({ 0.0, 0.0 }); });
    gradPhi.emplace_back([&](const Point& r)->Point { return Point({ sqrt(12.0 / hy / pow(hx, 3)), 0.0 }); });
    gradPhi.emplace_back([&](const Point& r)->Point { return Point({ 0.0, sqrt(12.0 / hx / pow(hy, 3)) }); });

    if (nShapes == 6)
    {
        phi.emplace_back([&](const Point& r)
            { return sqrt(180.0 / pow(hx, 5) / hy) * ( sqr(r.x() - center.x()) - sqr(hx) / 12.0); });
        phi.emplace_back([&](const Point& r)
            { return sqrt(144.0 / pow(hx, 3) / pow(hy, 3)) * (r.x() - center.x()) * (r.y() - center.y()); });
        phi.emplace_back([&](const Point& r)
            { return sqrt(180.0 / pow(hy, 5) / hx) * ( sqr(r.y() - center.y()) - sqr(hy) / 12.0); });

        gradPhi.emplace_back([&](const Point& r)->Point
            { return Point({ sqrt(180.0 / pow(hx, 5) / hy) * 2 * (r.x() - center.x()) , 0.0}); });
        gradPhi.emplace_back([&](const Point& r)->Point
            { return Point({ sqrt(144.0 / pow(hx, 3) / pow(hy, 3)) * (r.y() - center.y()), \
                             sqrt(144.0 / pow(hx, 3) / pow(hy, 3)) * (r.x() - center.x()) }); });
        gradPhi.emplace_back([&](const Point& r)->Point
            { return Point({ 0.0, sqrt(180.0 / pow(hy, 5) / hx) * 2.0 * (r.y() - center.y()) }); });

    }
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

vector<shared_ptr<Cell>> Cell::findNeighbourCellsX() const
{
    vector<shared_ptr<Cell>> neibCells;

    // for left edge
    if (edges[3]->neibCells.size() == 2)
        neibCells.emplace_back(edges[3]->neibCells[0]);

    // for right edge
    if (edges[1]->neibCells.size() == 2)
        neibCells.emplace_back(edges[1]->neibCells[1]);

    return neibCells;
}

vector<shared_ptr<Cell>> Cell::findNeighbourCellsY() const
{
    vector<shared_ptr<Cell>> neibCells;

    // for bottom edge
    if (edges[0]->neibCells.size() == 2)
        neibCells.emplace_back(edges[0]->neibCells[0]);

    // for top edge
    if (edges[2]->neibCells.size() == 2)
        neibCells.emplace_back(edges[2]->neibCells[1]);

    return neibCells;
}


// ----- RKDG -----

numvector<double, 5> Cell::reconstructSolution(const Point& point ) const
{
//    if (!insideCell(point))
//    {
//        std::cout << "Error: point (" << point.x() << ", " << point.y() << ") is not inside cell #" << number << std::endl;
//        std::cout << "Cell nodes:" << std::endl;

//        vector<shared_ptr<Point>> ccoord = getCellCoordinates();

//        for (size_t i = 0; i < ccoord.size(); ++i)
//            std::cout << "(" << ccoord[i]->x() << "; " << ccoord[i]->y() << ")" << endl;

//        exit(1);
//    }

    numvector<double, 5> sol(0.0);

    for (int i = 0; i < 5; ++i)
    {
        for (int j = 0; j < nShapes; ++j)
            sol[i] += phi[j](point) * problem->alpha[number][i * nShapes + j];
    }

    return sol;

} // end reconstructSolution


double Cell::reconstructSolution(const Point& point, int numSol ) const
{
//    if (!insideCell(point))
//    {
//        std::cout << "Error: point (" << point.x() << ", " << point.y() << ") is not inside cell #" << number << std::endl;
//        std::cout << "Cell nodes:" << std::endl;

//        vector<shared_ptr<Point>> ccoord = getCellCoordinates();

//        for (size_t i = 0; i < ccoord.size(); ++i)
//            std::cout << "(" << ccoord[i]->x() << "; " << ccoord[i]->y() << ")" << endl;

//        exit(1);
//    }

    double sol(0.0);
    
    for (int j = 0; j < nShapes; ++j)
        sol += phi[j](point) * problem->alpha[number][numSol * nShapes + j];

    return sol;
} // end reconstructSolution

numvector<double, 5 * nShapes> Cell::projection(std::function<numvector<double,5>(const Point& point)>& foo) const
{
    numvector<double, 5 * nShapes> alpha;

    for (int q = 0; q < nShapes; ++q)
    {
        std::function<numvector<double, 5>(const Point&)> f = \
                [&](const Point& p) {  return phi[q](p) * foo(p); };

        numvector<double, 5> buffer = integrate(f);

        for (int p = 0; p < 5; ++p)
        {
            alpha[p*nShapes + q] = buffer[p];
        }// for p
    }

    return alpha;

} // end projection

numvector<double, 5> Cell::integrate( const std::function<numvector<double, 5>(const Point &)>& f) const
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
} // end of cell integral

double Cell::getNormQ(int numSol) const
{
    vector<double> rhoGP(nGP);
    
    for (int i = 0; i < nGP; ++i)
	rhoGP[i] = reconstructSolution(gPoints2D[i],numSol);
    
    return *max_element(rhoGP.begin(), rhoGP.end());
}

