#pragma once

#include "numvector.h"
#include <vector>
#include <fstream>

namespace std
{

class Mesh2D
{
private:

    ofstream writer;

    

public:
	//- Space step in x direction
	double hx;

	//- Space step in y direction
	double hy;

    //- Number of internal cells
    double nInternalCells;

    //- Number of ghost cells
    double nGhostCells;



    //- Coordinates of nodes (x,y)
    vector<numvector<double, 2>> nodes;

    //- Horizontal edges (nodeLeft, nodeRight)
    vector<numvector<int, 2>> edgesHor;

    //- Vertical edges (nodeDown, nodeUp)
    vector<numvector<int, 2>> edgesVer;

    //- Mesh cells (edgeDown, edgeUp, edgeLeft, edgeRight)
    vector<numvector<int, 4>> cells;

    //- Internal mesh cells
    vector<numvector<int, 4>> internalCells;

    //- Ghost mesh cells (for boundary)
    vector<numvector<int, 4>> ghostCells;

	//- Cell centers
	vector<numvector<double, 2>> cellCenters;

	//- Neighbour cells for horizontal edges
	vector<numvector<int, 2>> neighbCellsHorEdges;

	//- Neighbour cells for vertical edges
	vector<numvector<int, 2>> neighbCellsVerEdges;

public:

    //- Default constructor
    Mesh2D() {}

    //- Copy constructor
    Mesh2D(const Mesh2D &mesh);

    //- operator =
    Mesh2D& operator=(const Mesh2D&);

    //- construct mesh by number of cells and size of flow domain
	Mesh2D(int nx, int ny, double Lx, double Ly);

    ~Mesh2D();

	//- get coordinates of nodes of given mesh cell
	numvector<numvector<double,2>,4> getCellCoordinates(int iCell);

    //- transform global coordinate system for give cell to local CS
    numvector<double, 2> globalToLocal(int iCell, numvector<double, 2> coord);

    //- transform global coordinate system for give cell to local CS
    numvector<double, 2> localToGlobal(int iCell, numvector<double, 2> coord);

    //- export mesh
    void exportMesh();

};// end Mesh 2D

} // end namespace std;

