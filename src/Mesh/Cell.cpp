#include "Cell.h"
#include "Edge.h"
#include <iostream>

using namespace std;

// ------------------ Constructors & Destructor ----------------

Cell::Cell(const std::vector<std::shared_ptr<Point> > &defNodes, const vector<std::shared_ptr<Edge>> &defEdges, const Problem& prb) : problem(prb)
{
    nodes = defNodes;
    edges = defEdges;

    nEntities = defNodes.size();
    
    //step.x() = edges[0]->nodes[1]->x() - edges[0]->nodes[0]->x();
    //step.y() = edges[nEdges-1]->nodes[1]->y() - edges[nEdges-1]->nodes[0]->y();
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

Point Cell::localToGlobal(const Point& localPoint) const
{
    vector<double> ff = {0,0,0};

    if (nEntities == 3)
        ff = {
                1.0 - localPoint.x() - localPoint.y(), \
                localPoint.x(), \
                localPoint.y() \
            };
    else if (nEntities == 4)
        ff = {
                0.25 * (1.0 - localPoint.x()) * (1.0 - localPoint.y()), \
                0.25 * (1.0 + localPoint.x()) * (1.0 - localPoint.y()), \
                0.25 * (1.0 + localPoint.x()) * (1.0 + localPoint.y()), \
                0.25 * (1.0 - localPoint.x()) * (1.0 + localPoint.y()) \
        };
    else
    {
        cout << "localToGlobal: Polygonal elements when nNodes > 4 are not supported\n";
        exit(0);
    }

    Point globalPoint({ 0.0, 0.0 });

    for (int i = 0; i < nEntities; ++i)
    {
        globalPoint[0] += nodes[i]->x() * ff[i];
        globalPoint[1] += nodes[i]->y() * ff[i];
    }

    return globalPoint;
}

void Cell::setArea()
{
    int n = nodes.size();

    area = 0.0;

    for (int i = 0; i < n-1; ++i)
        area += nodes[i]->x() * nodes[i+1]->y() - nodes[i]->y() * nodes[i+1]->x();

    area += nodes[n-1]->x() * nodes[0]->y() - nodes[n-1]->y() * nodes[0]->x();

    area = 0.5 * fabs(area);
}

void Cell::setJacobian()
{
    function<double(const Point&)> fJ;

    if (nEntities == 3)
    {
        fJ = [&](const Point& r){ return area; };
    }
    else if (nEntities == 4)
    {
        fJ = [&](const Point& r)
        {
            double dpde = (nodes[1]->x() - nodes[0]->x()) * (1.0 - r.y()) + (nodes[2]->x() - nodes[3]->x()) * (1.0 + r.y());
            double dqde = (nodes[1]->y() - nodes[0]->y()) * (1.0 - r.y()) + (nodes[2]->y() - nodes[3]->y()) * (1.0 + r.y());
            double dpdn = (nodes[3]->x() - nodes[0]->x()) * (1.0 - r.x()) + (nodes[2]->x() - nodes[1]->x()) * (1.0 + r.x());
            double dqdn = (nodes[3]->y() - nodes[0]->y()) * (1.0 - r.x()) + (nodes[2]->y() - nodes[1]->y()) * (1.0 + r.x());

            return 0.25 * (dpde * dqdn - dpdn * dqde);
        };
    }
    else
    {
        cout << "setJacobian: Polygonal elements when nNodes > 4 are not supported\n";
        exit(0);
    }

    J.resize(nGP);

    for (int i = 0; i < nGP; ++i)
        J[i] = fJ(gPoints2D[i]);
}

void Cell::setGaussPoints()
{    
    if (nShapes == 6)
    {
        // CHECK FOR UNSTRUCTURED MESHES!!!!

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
        if (nEntities == 4)
        {
            nGP = 4;

            double isqrt3 = 0.57735026918962576;

            for (int i = -1; i <= 1; i += 2)
                for (int j = -1; j <= 1; j += 2)
                    gPoints2D.push_back(localToGlobal(Point({ i * isqrt3, j * isqrt3 })));

            gWeights2D = { 1.0, 1.0, 1.0, 1.0 };
        }
        else if (nEntities == 3)
        {
//            nGP = 1;

//            gPoints2D.push_back(localToGlobal(Point({ 1.0/3.0, 1.0/3.0 })));

//            gWeights2D = { 1.0 };
            nGP = 3;

            gPoints2D.push_back(localToGlobal(Point({ 1.0/6.0, 1.0/6.0 })));
            gPoints2D.push_back(localToGlobal(Point({ 2.0/3.0, 1.0/6.0 })));
            gPoints2D.push_back(localToGlobal(Point({ 1.0/6.0, 2.0/3.0 })));

            gWeights2D = { 1.0/3.0, 1.0/3.0, 1.0/3.0 };
        }
        else
        {
            cout << "Cell #" << number << ", ";
            cout << "setGP: nNodes = " << nEntities << ", nNodes > 4 not supported\n";
            exit(0);
        }
    }
    else
    {
        cout << "Wrong number of basis functions " << nShapes;
        cout << "Available: 1,3,6 FF in 2D case \n";

        exit(0);
    }
}

void Cell::setBasisFunctions()
{
//    double isqrta  = 1.0 / sqrt(area);
//    double isqrta2 = 1.0 / sqrt(area * area / 6.0 / sqrt(3));

//    double isqrta  = 1.0;
//    double isqrta2 = 1.0;

    offsetPhi.push_back(1.0 / sqrt(area));
    offsetPhi.push_back(1.0 / sqrt(area * area / 6.0 / sqrt(3)));
    offsetPhi.push_back(1.0 / sqrt(area * area / 6.0 / sqrt(3)));

    phi.reserve(nShapes);
    gradPhi.reserve(nShapes);

//    phi.emplace_back([&](const Point& r){ return 1.0 ; });
//    phi.emplace_back([&](const Point& r){ return (r.x() - center.x()); });
//    phi.emplace_back([&](const Point& r){ return (r.y() - center.y()); });


    phi.emplace_back([&](const Point& r){ return 1.0 / sqrt(area); });
    phi.emplace_back([&](const Point& r){ return (r.x() - center.x()) * 1.0 / sqrt(area * area / 6.0 / sqrt(3)); });
    phi.emplace_back([&](const Point& r){ return (r.y() - center.y()) * 1.0 / sqrt(area * area / 6.0 / sqrt(3)); });

    gradPhi.emplace_back([&](const Point& r)->Point { return Point({ 0.0    , 0.0 }); });
    gradPhi.emplace_back([&](const Point& r)->Point { return Point({ 1.0 / sqrt(area * area / 6.0 / sqrt(3)), 0.0 }); });
    gradPhi.emplace_back([&](const Point& r)->Point { return Point({ 0.0    , 1.0 / sqrt(area * area / 6.0 / sqrt(3)) }); });

    if (nShapes == 6)
    {
        // CHECK FOR UNSTRUCTURED MESHES!!!
        //just for compilation, not used!!!
        double hx = 1.0;
        double hy = 1.0;



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

void Cell::setNonOrthoMatrix()
{
    nonOrthoMatrix.resize(nShapes);

    for (int i = 0; i < nShapes; ++i)
        nonOrthoMatrix[i].resize(nShapes);

    std::function<double(const Point&)> f ;

    for (int i = 0; i < nShapes; ++i)
    {
        f = [&](const Point& p) { return phi[i](p) * phi[i](p); };
        nonOrthoMatrix[i][i] = integrate(f);

        for (int j = i+1; j < nShapes; ++j)
        {
            f = [&](const Point& p) {  return phi[i](p) * phi[j](p); };
            nonOrthoMatrix[i][j] = integrate(f);

            nonOrthoMatrix[j][i] = nonOrthoMatrix[i][j];
        }
    }
}

// ------------------ Public class methods ---------------------

// ----- geometric -----

vector<shared_ptr<Point>> Cell::getCellCoordinates() const
{
    vector<shared_ptr<Point>> nodeCoordinates(nEntities);

    for (int i = 0; i < nEntities; ++i)
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
    numvector<double, 5> sol(0.0);

    for (int i = 0; i < 5; ++i)
    {
        for (int j = 0; j < nShapes; ++j)
            sol[i] += phi[j](point) * problem.alpha[number][i * nShapes + j];
    }

    return sol;

} // end reconstructSolution


double Cell::reconstructSolution(const Point& point, int numSol ) const
{
    double sol(0.0);
    
    for (int j = 0; j < nShapes; ++j)
        sol += phi[j](point) * problem.alpha[number][numSol * nShapes + j];

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

double Cell::integrate( const std::function<double(const Point &)>& f) const
{
    double res = 0.0;

    for (int i = 0; i < nGP; ++i)
    {
        double resF = f(gPoints2D[i]);

        res += gWeights2D[i] * resF * J[i];

    }

    return res;

} // end integrate 2D of vector function

numvector<double, 5> Cell::integrate( const std::function<numvector<double, 5>(const Point &)>& f) const
{
    numvector<double, 5> res = 0.0;

    for (int i = 0; i < nGP; ++i)
    {
         numvector<double,5> resF = f(gPoints2D[i]);

        for (int k = 0; k < 5; ++k)
        {
           res[k] += gWeights2D[i] * resF[k] * J[i];
        }
    }

    return res;

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
            resV = problem.fluxF(sol) * gradPhi[q](gPoints2D[i])[0] + \
                   problem.fluxG(sol) * gradPhi[q](gPoints2D[i])[1];

            for (int p = 0; p < 5; ++p)
                res[p * nShapes + q] += resV[p] * gWeights2D[i] * J[i];
        }
    }

    return res;
} // end of cell integral

double Cell::getNormQ(int numSol) const
{
    vector<double> rhoGP(nGP);
    
    for (int i = 0; i < nGP; ++i)
	rhoGP[i] = reconstructSolution(gPoints2D[i],numSol);
    
    return *max_element(rhoGP.begin(), rhoGP.end());
}

numvector<double, 5 * nShapes> Cell::correctNonOrtho(const numvector<double, 5 * nShapes>& rhs) const
{
    numvector<double, 5 * nShapes> alphaCorr;

//    vector<vector<double>> slae(nShapes);

//    for (int i = 0; i < nShapes; ++i)
//        slae[i].resize(nShapes+1);

    for (int iSol = 0; iSol < 5; ++iSol)
    {
        // prepare data for gaussian solver
//        for (int i = 0; i < nShapes; ++i)
//        {
//            for (int j = 0; j < nShapes; ++j)
//                slae[i][j] = nonOrthoMatrix[i][j];
//            slae[i][nShapes] = rhs[i + iSol*nShapes];
//        }

        // solve slae
        vector<double> solution(nShapes); //for 3 ff!!!

        double ref = nonOrthoMatrix[1][1] / nonOrthoMatrix[2][1];

        solution[0] = rhs[iSol*nShapes];
        solution[2] = (rhs[iSol*nShapes + 2] * ref - rhs[iSol*nShapes + 1]) \
                / (nonOrthoMatrix[2][2] * ref);
        solution[1] = (rhs[iSol*nShapes + 1] - solution[2] * nonOrthoMatrix[1][2]) \
                / (nonOrthoMatrix[1][1]);
        //set solution to appropriate positions
        for (int i = 0; i < nShapes; ++i)
            alphaCorr[i + iSol*nShapes] = solution[i];
    }

    return alphaCorr;
}

