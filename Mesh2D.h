#pragma once

#include "numvector.h"
#include "Edge.h"
#include "RKDGCell.h"
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
    int nInternalCells;

    //- Number of ghost cells
    int nGhostCells;



    //- Coordinates of nodes (x,y)
    vector<numvector<double, 2>> nodes;

    //- Horizontal edges (nodeLeft, nodeRight)
    vector<Edge*> edgesHor;

    //- Vertical edges (nodeDown, nodeUp)
    vector<Edge*> edgesVer;

    //- Mesh cells (edgeDown, edgeUp, edgeLeft, edgeRight)
    vector<RKDGCell*> cells;

    //- Internal mesh cells
    vector<RKDGCell*> internalCells;

    //- Ghost mesh cells (for boundary)
    vector<RKDGCell*> ghostCells;

public:

    //- Default constructor
    Mesh2D() {}

    //- Construct mesh by number of cells and size of flow domain
	Mesh2D(int nx, int ny, double Lx, double Ly);

    ~Mesh2D();

    //- Get coordinates of nodes of given mesh cell
	numvector<numvector<double,2>,4> getCellCoordinates(int iCell);

    //- Export mesh
    void exportMesh();

};// end Mesh 2D

} // end namespace std;

