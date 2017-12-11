#include "Mesh2D.h"
#include <iostream>

using namespace std;

// ------------------ Constructors & Destructors ----------------


// nx --- number of cells along x axis
// ny --- number of cells along y axis
Mesh2D::Mesh2D(int nx, int ny, double Lx, double Ly)
{
    // define number of nodes, edges, cells

    int nNodes = (nx + 1)*(ny + 1);
	int nEdgesHor = nx*(ny + 1);
	int nEdgesVer = ny*(nx + 1);

    nCells = nx*ny;
    //nGhostCells = 2 * nx + 2 * ny;

    // space step in x, y direction

    hx = Lx / nx;
    hy = Ly / ny;


    // reserve memory
	
	nodes.reserve(nNodes);
    edgesHor.reserve(nEdgesHor);
    edgesVer.reserve(nEdgesVer);
//    cells.reserve(nInternalCells);
    //ghostCells.reserve(nGhostCells);
    //cells.reserve(nInternalCells+nGhostCells);

    // fill nodes

	for (int i = 0; i < ny + 1; ++i)
    {
		for (int j = 0; j < nx + 1; ++j)
        {
            Point node(j*hx, i*hy);
            nodes.push_back(node);
        }
    }

//    //for (int i = 0; i < nNodes; ++i)
//    //    cout << nodes[i].x() << ' ' << nodes[i].y() << endl;


    // get horizontal edges

    for (int i = 0; i < ny + 1; ++i)
    {
        EdgeBoundaryInfty currentEdgeLeft (&nodes[i*(nx + 1)], &nodes[i*(nx + 1) + 1]);
        edgesHor.push_back(currentEdgeLeft);
        edgesBound.push_back(currentEdgeLeft);

        for (int j = 1; j < nx - 1; ++j)
        {
            EdgeInternal currentEdge (&nodes[i*(nx + 1) + j], &nodes[i*(nx + 1) + j + 1]);
            edgesHor.push_back(currentEdge);
        }

        EdgeBoundaryInfty currentEdgeRight (&nodes[i*(nx + 1) + nx - 1], &nodes[i*(nx + 1) + nx]);
        edgesHor.push_back(currentEdgeRight);
        edgesBound.push_back(currentEdgeRight);
    }

    // get vertical edges

    for (int i = 0; i < ny; ++i)
    {
        EdgeBoundaryInfty currentEdgeLeft (&nodes[i*(nx + 1)], &nodes[i*(nx + 1) + nx + 1]);
        edgesVer.push_back(currentEdgeLeft);
        edgesBound.push_back(currentEdgeLeft);

        for (int j = 0; j < nx + 1; ++j)
        {
            EdgeInternal currentEdge (&nodes[i*(nx + 1) + j], &nodes[i*(nx + 1) + j + nx + 1]);
            edgesVer.push_back(currentEdge);
        }

        EdgeBoundaryInfty currentEdgeRight (&nodes[i*(nx + 1) + nx - 1], &nodes[i*(nx + 1) + nx + nx]);
        edgesVer.push_back(currentEdgeRight);
        edgesBound.push_back(currentEdgeRight);
    }

    // get cells as edges (counter-clockwise: lb -> rb -> ru -> rl)
    for (int i = 0; i < ny; ++i)
    {
        for (int j = 0; j < nx; ++j)
        {
            //define list of edges

            numvector<Edge*,4> edges = {&edgesHor[i*(nx)+j], \
                                        &edgesVer[(i)*(nx + 1) + j + 1], \
                                        &edgesHor[(i + 1)*(nx)+j], \
                                        &edgesVer[(i)*(nx + 1) + j] };

            Cell currentCell(edges);


            // this cell is neighbour for its edges

            edgesHor[i       * nx       + j].neibCells.push_back(&currentCell);
            edgesVer[i       * (nx + 1) + j + 1].neibCells.push_back(&currentCell);
            edgesHor[(i + 1) * nx       + j].neibCells.push_back(&currentCell);
            edgesVer[i       * (nx + 1) + j].neibCells.push_back(&currentCell);

            currentCell.number = i;

            // add cells in list

            cells.push_back(currentCell);
        }
    }

    writer.open("mesh2D");
}

Mesh2D::~Mesh2D()
{
    writer.close();
}

// ------------------ Private class methods --------------------

// ------------------ Public class methods ---------------------


void Mesh2D::exportMesh()
{
    writer << "hx = " << hx << endl;
    writer << "hy = " << hy << endl;

    writer << "cellCenters" << endl;

    for (int i = 0; i < nCells; ++i)
    {
        writer << cells[i].center.x() << ' ' << cells[i].center.y() << endl;
    }

    cout << "Mesh export OK" << endl;
}
