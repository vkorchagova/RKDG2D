#include "Mesh2D.h"
#include <iostream>

using namespace std;

// ------------------ Constructors & Destructors ----------------

Mesh2D::Mesh2D(const Mesh2D &mesh)
{
    hx = mesh.hx;
    hy = mesh.hy;
    nodes = mesh.nodes;
    edgesHor = mesh.edgesHor;
    edgesVer = mesh.edgesVer;
    cells = mesh.cells;
    cellCenters = mesh.cellCenters;
    neighbCellsHorEdges = mesh.neighbCellsHorEdges;
    neighbCellsVerEdges = mesh.neighbCellsVerEdges;
}

// nx --- number of cells along x axis
// ny --- number of cells along y axis
Mesh2D::Mesh2D(int nx, int ny, double Lx, double Ly)
{
	int nNodes = (nx + 1)*(ny + 1);
	int nEdgesHor = nx*(ny + 1);
	int nEdgesVer = ny*(nx + 1);
    nInternalCells = nx*ny;
    nGhostCells = 2 * nx + 2 * ny;

	hx = Lx / nx;
	hy = Ly / ny;
	
	nodes.reserve(nNodes);
	edgesHor.reserve(nEdgesHor);
	edgesVer.reserve(nEdgesVer);
    internalCells.reserve(nInternalCells);
    ghostCells.reserve(nGhostCells);
    cells.reserve(nInternalCells+nGhostCells);
    cellCenters.reserve(nInternalCells+nGhostCells);

	// fill nodes
	for (int i = 0; i < ny + 1; ++i)
		for (int j = 0; j < nx + 1; ++j)
			nodes.push_back({ j*hx, i*hy });


	//get edges
	for (int i = 0; i < ny + 1; ++i)
		for (int j = 0; j < nx; ++j)
			edgesHor.push_back({ i*(nx + 1) + j, i*(nx + 1) + j + 1 });


	for (int i = 0; i < ny; ++i)
		for (int j = 0; j < nx + 1; ++j)
			edgesVer.push_back({ i*(nx + 1) + j, i*(nx + 1) + j + nx + 1 });

	// get cells as edges
    for (int i = 0; i < ny; ++i)
        for (int j = 0; j < nx; ++j)
            internalCells.push_back({i*(nx)+j, (i + 1)*(nx)+j, (i)*(nx + 1) + j, (i)*(nx + 1) + j + 1 });

	// get ghost cells

	//down
	for (int i = 0; i < nx; ++i)
        ghostCells.push_back({ -1, i, -1, -1 });

	//up
	for (int i = 0; i < nx; ++i)
        ghostCells.push_back({ nx*ny + i, -1, -1, -1 });

	//left
	for (int j = 0; j <  ny; ++j)
        ghostCells.push_back({ -1, -1, -1, (j)*(nx + 1) });
	 
	//right
	for (int j = 0 ; j < ny; ++j)
        ghostCells.push_back({ -1, -1, (j)*(nx + 1) + nx, -1 });

    cells.insert(cells.end(),internalCells.begin(),internalCells.end());
    cells.insert(cells.end(),ghostCells.begin(),ghostCells.end());

	//get cell centers
	numvector<double, 2> offset = {0.5*hx, 0.5*hy};

	//for internal cells
    for (int i = 0; i < nInternalCells; ++i)
	{
		numvector<double,2> lowerLeftNode = getCellCoordinates(i)[0];
		cellCenters.push_back(lowerLeftNode + offset);
	}

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

	// get neighbour cells for horizontal edges
	for (int i = 0; i < nEdgesHor; ++i)
	{
        numvector<int, 2> neighbours = {-1,-1};

        for (int j = 0; j < nInternalCells + nGhostCells; ++j)
		{

            if (i == cells[j][0])
				neighbours[1] = j;
			if (i == cells[j][1])
				neighbours[0] = j;
			if (neighbours[0] > -1 && neighbours[1] > -1)
				break;
		}

		neighbCellsHorEdges.push_back(neighbours);
	}

	// get neighbour cells for vertical edges
	for (int i = 0; i < nEdgesVer; ++i)
	{
		numvector<int, 2> neighbours = { -1, -1 };

        // more optimal search?..
        for (int j = 0; j < nInternalCells + nGhostCells; ++j)
		{
			if (i == cells[j][2])
				neighbours[1] = j;
			if (i == cells[j][3])
				neighbours[0] = j;
			if (neighbours[0] > -1 && neighbours[1] > -1)
				break;
		}

		neighbCellsVerEdges.push_back(neighbours);
	}
}

// ------------------ Private class methods --------------------

// ------------------ Public class methods ---------------------

numvector<numvector<double, 2>, 4> Mesh2D::getCellCoordinates(int iCell)
{
	int iEdgeDown = cells[iCell][0];
	int iEdgeUp = cells[iCell][1];

	numvector<int, 2> iNodesDown = edgesHor[iEdgeDown];
	numvector<int, 2> iNodesUp = edgesHor[iEdgeUp];

	numvector<numvector<double, 2>, 4> nodeCoordinates;

	nodeCoordinates[0] = nodes[iNodesDown[0]];
	nodeCoordinates[1] = nodes[iNodesDown[1]];
	nodeCoordinates[2] = nodes[iNodesUp[1]];
	nodeCoordinates[3] = nodes[iNodesUp[0]];

	return nodeCoordinates;
}

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
