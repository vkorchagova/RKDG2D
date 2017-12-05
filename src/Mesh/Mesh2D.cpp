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
    internalCells.reserve(nInternalCells);
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

            internalCells.push_back(currentCell);
        }
    }

    /*
    for (int i = 0; i < nInternalCells; ++i)
        for (int j = 0; j < 4; ++j)
        {
            cout << "edge num = " << j << endl;
            cout << internalCells[i].edges[j]->nodes[0]->x() << ' ' \
                 << internalCells[i].edges[j]->nodes[0]->y() << " to " \
                 << internalCells[i].edges[j]->nodes[1]->x() << ' ' \
                 << internalCells[i].edges[j]->nodes[1]->y() << endl;
        }

        */

    // get ghost cells
    // not needed? Boundary conditions will be used
    /*
	//down
	for (int i = 0; i < nx; ++i)
    {
        numvector<Edge*,4> edges = {NULL, \
                                    NULL, \
                                    NULL, \
                                    edgesHor[i] };
        Cell currentCell(edges);

        ghostCells.push_back(&currentCell);
    }

	//up
	for (int i = 0; i < nx; ++i)
    {
        numvector<Edge*,4> edges = {edgesHor[nx*ny + i], \
                                    NULL, \
                                    NULL, \
                                    NULL };

        Cell currentCell(edges);

        ghostCells.push_back(&currentCell);
    }

	//left
	for (int j = 0; j <  ny; ++j)
    {
        numvector<Edge*,4> edges = {NULL, \
                                    NULL, \
                                    edgesVer[j*(nx + 1)], \
                                    NULL};

        Cell currentCell(edges);

        ghostCells.push_back(&currentCell);
    }
	 
	//right
	for (int j = 0 ; j < ny; ++j)
    {
        numvector<Edge*,4> edges = {NULL, \
                                    edgesVer[j*(nx + 1) + nx], \
                                    NULL, \
                                    NULL};

        Cell currentCell(edges);

        ghostCells.push_back(&currentCell);
    }

*/
    cells.insert(cells.end(),internalCells.begin(),internalCells.end());
   //cells.insert(cells.end(),ghostCells.begin(),ghostCells.end());


    /*

	//for ghost cells
	for (int i = 0; i < nx; ++i)
	{
		numvector<double, 2> currentCellCenter = {(0.5+i)*hx, -0.5*hy};
		cellCenters.push_back(currentCellCenter);
	}

	for (int i = 0; i < nx; ++i)
	{
		numvector<double, 2> currentCellCenter = { (0.5 + i)*hx, Ly + 0.5*hy };
		cellCenters.push_back(currentCellCenter);
	}

	for (int j = 0; j < ny; ++j)
	{
		numvector<double, 2> currentCellCenter = { -0.5*hx, (j + 0.5)*hy };
		cellCenters.push_back(currentCellCenter);
	}

	for (int j = 0; j < ny; ++j)
	{
		numvector<double, 2> currentCellCenter = { Lx + 0.5*hx, (j + 0.5)*hy };
		cellCenters.push_back(currentCellCenter);
	}
    */

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


/*
numvector<double, 2> Mesh2D::globalToLocal(int iCell, numvector<double, 2> coord)
{
	//numvector<double, 2> steps = {1.0/hx, 1.0/hy};
	//return 2 * (coord - cellCenters[iCell]) * steps;

	numvector<double, 2> offset = coord - cellCenters[iCell];
	
	return{ 2 * offset[0] / hx, 2 * offset[1] / hy };
}

numvector<double, 2> Mesh2D::localToGlobal(int iCell, numvector<double, 2> coord)
{
	numvector<double, 2> steps = { 0.5*coord[0] * hx, 0.5*coord[1] * hy };
	return cellCenters[iCell] + steps;
}
*/

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
