#ifndef CELL_H
#define CELL_H

/// ------------------------------
/// Cell class
/// ------------------------------
/// Consists of 4 edges
/// Know its Gauss points
/// Can evaluate RHS inside cell
/// ------------------------------

#include "Point.h"
//#include "Edge.h"
#include "numvector.h"
//#include "Problem.h"
#include "defs.h"

#include <functional>
#include <memory>
#include <math.h>
#include <algorithm>

class Edge;
class Problem;


class Cell
{

private:

    //- Number of Gauss points
    int nGP;

    //- Gauss points
    std::vector<Point> gPoints2D;

    //- Gauss weights
    std::vector<double> gWeights2D;

    //- Jacobian
    std::vector<double> J;

    //- Area of cell
    double area;

    //- Mass center of cell
    Point center;


    //- Gramian matrix
    std::vector<std::vector<double>> gramian;

    //- local [-1,1]x[-1,1] to global rectangular cell
    Point localToGlobal(const Point& localPoint) const;




public:

    /// geometric variables

    //- Number of cell
    int number;

    //- Number of entities (edges or numbers - mo matter)
    int nEntities;

    //- Nodes in cell
    std::vector<std::shared_ptr<Node>> nodes;

    //- Edges in cell
    std::vector<std::shared_ptr<Edge>> edges;

    //- Neighbour cells
    std::vector<std::shared_ptr<Cell>> neibCells;

    //- Return area of cell
    double getArea() const { return area; }

    //- Return center of cell
    const Point& getCellCenter() const { return center; }

    //- Calculate element area
    void setArea();

    //- Set cell center
    void setCellCenter();

    //- Define Jacobian function
    void setJacobian();

    //- Set Gauss points
    void setGaussPoints();

    //- Set basis functions
    void setBasisFunctions();

    /// RKDG variables

    //- Pointer to problem
    const Problem& problem;

    //- List of basis functions
    std::vector<std::function<double(const Point&)>> phi;

    //- Offset for basis functions
    //debug
    std::vector<double> offsetPhi;

    //- Gradient of basis functions
    std::vector<std::function<Point(const Point&)>> gradPhi;

public:

    Cell(const Problem& prb) : problem(prb) {}

    //- Construct cell using vectors of nodes and edges
    Cell(const std::vector<std::shared_ptr<Node> > &nodes, const std::vector<std::shared_ptr<Edge> > &edges, const Problem& prb);

    //- Destructor
    ~Cell() {}

    /// geometric methods

    //- Calculate coordinates of cell nodes
    std::vector<std::shared_ptr<Point>> getCellCoordinates() const;

    //- Find neighbour cells
    void findNeighbourCells();

    //- Find neighbour cells in X direction
    //DELETE
    std::vector<std::shared_ptr<Cell>> findNeighbourCellsX() const;

    //- Find neighbour cells in Y direction
    //DELETE
    std::vector<std::shared_ptr<Cell>> findNeighbourCellsY() const;

    //- Check if point belongs cell
    bool insideCell(const Point& point) const;

    //- TODO: list of neighbour cells


    /// RKDG methods

    //- Reconstruct solution (get coeffs from problem)
    numvector<double, 5> reconstructSolution(const Point& point) const;
    double reconstructSolution(const Point& point, int numSol) const;

    numvector<double, 5> reconstructSolution(const std::shared_ptr<Point> point) const
        { return reconstructSolution(*point); }
    double reconstructSolution(const std::shared_ptr<Point> point, int numSol) const
        { return reconstructSolution(*point, numSol); }

    //- Reconstruct solution using given coeffs
    numvector<double, 5> reconstructSolution(const Point& point, const numvector<double, 5*nShapes>& alpha) const;
    double reconstructSolution(const Point& point, const numvector<double, 5*nShapes>& alpha, int numSol) const;

    numvector<double, 5> reconstructSolution(const std::shared_ptr<Point> point, const numvector<double, 5*nShapes>& alpha) const
        { return reconstructSolution(*point, alpha); }
    double reconstructSolution(const std::shared_ptr<Point> point, const numvector<double, 5*nShapes>& alpha, int numSol) const
        { return reconstructSolution(*point, alpha, numSol); }

    //- Reconstruct coefficients using Riemann invariants
    numvector<double, 5 * nShapes> reconstructCoefficients(const std::pair<numvector<double, 5 * nShapes>, numvector<double, 5 * nShapes>>& rI) const;
    numvector<double, 5 * nShapes> reconstructCoefficients(const numvector<double, 5 * nShapes>& rI, const Point& n) const;

    //- Get coefficients of projection of function foo onto cell basis
    numvector<double, 5 * nShapes> projection(std::function<numvector<double,5>(const Point& point)>& init) const;

    //- Get Riemann invariants using actual solution on cell
    std::pair<numvector<double, 5 * nShapes>, numvector<double, 5 * nShapes>> getRiemannInvariants();
    numvector<double, 5 * nShapes> getRiemannInvariants(const Point& n);

    //- Calculate Gramian matrix
    void setGramian();

    //- Solve SLAE in case of non-orthogonal functions
    numvector<double, 5 * nShapes> correctNonOrtho(const numvector<double, 5 * nShapes>& rhs) const;

    //- Reconstruct SLAE RHS after limitation in case of non-orthogonal functions
    numvector<double, 5 * nShapes> correctPrevIter(const numvector<double, 5 * nShapes>& rhs) const;

    //- Calculate \int_{cell} F(U) \nabla \phi_x + G(U) \nabla \phi_y
    numvector<double, 5 * nShapes> cellIntegral();

    //- Calculate total mass flux
    double totalMassFlux() const;

    //- Calculate total mass in cell
    double totalMass() const;

    //- 2D Gauss integration of scalar function
    double integrate( const std::function<double(const Point &)>& f) const;


    //- 2D Gauss integration of vector function
    numvector<double,5> integrate( const std::function<numvector<double, 5>(const Point&)>& f) const;
    
    //- Get norm of solution for indicators
    double getNormQ(int numSol) const;

    //- Get norm of pressure for indicators
    double getNormQp () const;
};

#endif // CELL_H
