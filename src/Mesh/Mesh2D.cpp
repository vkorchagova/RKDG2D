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

    nInternalCells = nx*ny;
    //nGhostCells = 2 * nx + 2 * ny;

    // space step in x, y direction

    hx = Lx / nx;
    hy = Ly / ny;


    // reserve memory
	
	nodes.reserve(nNodes);
	edgesHor.reserve(nEdgesHor);
	edgesVer.reserve(nEdgesVer);
    cells.reserve(nInternalCells);
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

    //for (int i = 0; i < nNodes; ++i)
    //    cout << nodes[i].x() << ' ' << nodes[i].y() << endl;


    // get horizontal edges

	for (int i = 0; i < ny + 1; ++i)
    {
		for (int j = 0; j < nx; ++j)
        {
            Edge currentEdge;
            currentEdge.nodes[0] = &nodes[i*(nx + 1) + j];
            currentEdge.nodes[1] = &nodes[i*(nx + 1) + j + 1];
            edgesHor.push_back(currentEdge);
        }
    }

    // get vertical edges

	for (int i = 0; i < ny; ++i)
    {
		for (int j = 0; j < nx + 1; ++j)
        {
            Edge currentEdge;
            currentEdge.nodes[0] = &nodes[i*(nx + 1) + j];
            currentEdge.nodes[1] = &nodes[i*(nx + 1) + j + nx + 1];
            edgesVer.push_back(currentEdge);
        }
    }

    // get cells as edges (counter-clockwise: lb -> rb -> ru -> rl)
    for (int i = 0; i < ny; ++i)
    {
        for (int j = 0; j < nx; ++j)
        {
            numvector<Edge*,4> edges = {&edgesHor[i*(nx)+j], \
                                        &edgesVer[(i)*(nx + 1) + j + 1], \
                                        &edgesHor[(i + 1)*(nx)+j], \
                                        &edgesVer[(i)*(nx + 1) + j] };

            Cell currentCell(edges);

            cells.push_back(currentCell);
        }
    }

    // get neighbour cells


    for (int i = 0; i < nInternalCells; ++i)
    {
        cells[i].number = i;

        //down
        cells[i].edges[0]->neibCells[1] = &cells[i];

        //right
        cells[i].edges[1]->neibCells[1] = &cells[i];

        //up
        cells[i].edges[2]->neibCells[0] = &cells[i];

        //left
        cells[i].edges[3]->neibCells[0] = &cells[i];
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

    for (int i = 0; i < nInternalCells; ++i)
    {
        writer << cells[i].center.x() << ' ' << cells[i].center.y() << endl;
    }

    cout << "Mesh export OK" << endl;
}
