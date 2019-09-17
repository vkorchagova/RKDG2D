#include "Cell.h"
#include "defs.h"
//#include "Edge.h"
#include "compService.h"
#include <iostream>

using namespace std;

// ------------------ Constructors & Destructor ----------------

//Cell::Cell(const std::vector<std::shared_ptr<Point>> &defNodes, const vector<std::shared_ptr<Edge>> &defEdges)
Cell::Cell(const std::vector<std::shared_ptr<Point>> &defNodes, const vector<std::shared_ptr<Edge>> &defEdges)
: nodes(defNodes), edges(defEdges)
{
    nEntities = defNodes.size();
    setArea();
    setGaussPoints();
    setJacobian();
    setCellCenter();
}

// ------------------ Private class methods --------------------


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
    vector<Point> localGP;

//    function<double(double)> chop = [](double num) { if (fabs(num) < 1e-12) return 0.0; return num;};

    if (nEntities == 3)
    {
        fJ = [&](const Point& r){ return area; };


        localGP.push_back(Point({ 1.0/6.0, 1.0/6.0 }));
        localGP.push_back(Point({ 2.0/3.0, 1.0/6.0 }));
        localGP.push_back(Point({ 1.0/6.0, 2.0/3.0 }));

    }
    else if (nEntities == 4)
    {
        fJ = [&](const Point& r)
        {
            double dpde = (nodes[1]->x() - nodes[0]->x()) * (1.0 - r.y()) + (nodes[2]->x() - nodes[3]->x()) * (1.0 + r.y());
            double dqde = (nodes[1]->y() - nodes[0]->y()) * (1.0 - r.y()) + (nodes[2]->y() - nodes[3]->y()) * (1.0 + r.y());
            double dpdn = (nodes[3]->x() - nodes[0]->x()) * (1.0 - r.x()) + (nodes[2]->x() - nodes[1]->x()) * (1.0 + r.x());
            double dqdn = (nodes[3]->y() - nodes[0]->y()) * (1.0 - r.x()) + (nodes[2]->y() - nodes[1]->y()) * (1.0 + r.x());

            return 0.0625 * fabs(dpde * dqdn - dpdn * dqde);
        };

        double isqrt3 = 0.57735026918962576;

        for (int i = -1; i <= 1; i += 2)
            for (int j = -1; j <= 1; j += 2)
                localGP.push_back(Point({ i * isqrt3, j * isqrt3 }));

    }
    else
    {
        cout << "setJacobian: Polygonal elements when nNodes > 4 are not supported\n";
        exit(0);
    }

    J.resize(nGP);

    for (int i = 0; i < nGP; ++i)
        J[i] = fJ(localGP[i]);
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

            // pure Gauss
            nGP = 4;

            double isqrt3 = 0.57735026918962576;

            for (int i = -1; i <= 1; i += 2)
                for (int j = -1; j <= 1; j += 2)
                   gPoints2D.push_back(localToGlobal(Point({ i * isqrt3, j * isqrt3 })));

            gWeights2D = { 1.0, 1.0, 1.0, 1.0 };

            // Gauss --- Lobatto

//            nGP = 9;

//            for (int i = -1; i <= 1; i++)
//                for (int j = -1; j <= 1; j++)
//                    gPoints2D.push_back(localToGlobal(Point({ i, j })));

//            const double frac43 = 1.3333333333333333;
//            const double frac13 = 0.3333333333333333;

//            vector <double> gWeights1D = { frac13, frac43, frac13 };

//            for (int i = 0; i <= 2; ++i)
//                for (int j = 0; j <= 2; ++j)
//                    gWeights2D.push_back(gWeights1D[i] * gWeights1D[j]);


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
            //cout << "Cell #" << number << ", ";
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

void Cell::setCellCenter()
{
    double xc = integrate(*this, function<double(const Point&)>([](const Point& r) {return r.x();}));
    double yc = integrate(*this, function<double(const Point&)>([](const Point& r) {return r.y();}));

    center.x() = xc / area;
    center.y() = yc / area;
}


// ------------------ Public class methods ---------------------

// ----- geometric -----

vector<Point> Cell::getCellCoordinates() const
{
    vector<Point> nodeCoordinates(nEntities);

    for (int i = 0; i < nEntities; ++i)
    {
        nodeCoordinates[i] = *(edges[i]->nodes[0]);
    }

    return nodeCoordinates;

} // end getCellCoordinates



bool Cell::inNeibVertList(const std::shared_ptr<Cell>& cell, std::vector<std::shared_ptr<Cell>>& neibVC)
{
    for (auto c : neibVC)
    {
        if (c->number == cell->number)
            return true;
    }
    return false;
}

bool Cell::hasCommonNode(const std::shared_ptr<Cell>& c1)
{
    for (auto n1 : c1->nodes)
        for (auto n2 : nodes)
            if (n1 == n2)
                return true;
    return false;
}
