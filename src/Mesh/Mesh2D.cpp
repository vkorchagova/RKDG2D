#include "Mesh2D.h"
#include <iostream>

using namespace std;

typedef EdgeBoundaryOpen edgeBoundaryT;

// ------------------ Constructors & Destructors ----------------


// nx --- number of cells along x axis
// ny --- number of cells along y axis
Mesh2D::Mesh2D(int nx, int ny, double Lx, double Ly) : nx(nx), ny(ny), Lx(Lx), Ly(Ly)
{
    // define number of nodes, edges, cells

    int nNodes = (nx + 1)*(ny + 1);
    int nEdgesHor = nx*(ny + 1);
    int nEdgesVer = ny*(nx + 1);

    nCells = nx*ny;

    // space step in x, y direction

    double hx = Lx / nx;
    double hy = Ly / ny;

    // reserve memory
	
    nodes.reserve(nNodes);
    edgesHor.reserve(nEdgesHor);
    edgesVer.reserve(nEdgesVer);
    
    // fill nodes

    for (int i = 0; i < ny + 1; ++i)
        for (int j = 0; j < nx + 1; ++j)
            nodes.emplace_back( Point({ j * hx, i * hy }) );


    // get horizontal edges

    for (int j = 0; j < nx; ++j)
        edgesHor.emplace_back( make_shared<edgeBoundaryT>(nodes[j], nodes[j + 1]) );

    for (int i = 1; i < ny; ++i)
       for (int j = 0; j < nx; ++j)
            edgesHor.emplace_back( make_shared<EdgeInternal>(nodes[i*(nx + 1) + j], nodes[i*(nx + 1) + j + 1]) );


    for (int j = 0; j < nx; ++j)
        edgesHor.emplace_back( make_shared<edgeBoundaryT>(nodes[ny*(nx + 1) + j], nodes[ny*(nx + 1) + j + 1]) );


    // get vertical edges

    for (int i = 0; i < ny; ++i)
    {
        edgesVer.emplace_back( make_shared<edgeBoundaryT>(nodes[i*(nx + 1)], nodes[i*(nx + 1) + nx + 1]) );

        for (int j = 1; j < nx; ++j)
            edgesVer.emplace_back( make_shared<EdgeInternal>(nodes[i*(nx + 1) + j], nodes[i*(nx + 1) + j + nx + 1]) );

        edgesVer.emplace_back( make_shared<edgeBoundaryT>(nodes[i*(nx + 1) + nx ], nodes[i*(nx + 1) + nx + nx + 1]) );
    }

    // set normals to horizontal edges
    for (int i = 0; i < nx; ++i)
        edgesHor[i]->n = Point({ 0.0, -1.0 });

    for (int i = nx; i < nEdgesHor; ++i)
        edgesHor[i]->n = Point({ 0.0,  1.0 });

    // set normals to vertical edges
    for (int i = 0; i < nEdgesVer; ++i)
        edgesVer[i]->n = (i % (nx + 1) == 0) ? Point({ -1.0, 0.0 }) : Point({ 1.0, 0.0 });


    // get cells as edges (counter-clockwise: lb -> rb -> ru -> rl)
    for (int i = 0; i < ny; ++i)
    {
        for (int j = 0; j < nx; ++j)
        {
            //define list of edges

            numvector<shared_ptr<Edge>,4> edges = {edgesHor[i*(nx)+j], \
                                                   edgesVer[(i)*(nx + 1) + j + 1], \
                                                   edgesHor[(i + 1)*(nx)+j], \
                                                   edgesVer[(i)*(nx + 1) + j] };

            // add cells in list
            cells.emplace_back( make_shared<Cell>(edges));

            // this cell is neighbour for its edges

            edgesHor[ i       *  nx      + j    ]->neibCells.push_back(cells[i*(nx)+j]);
            edgesVer[ i       * (nx + 1) + j + 1]->neibCells.push_back(cells[i*(nx)+j]);
            edgesHor[(i + 1)  *  nx      + j    ]->neibCells.push_back(cells[i*(nx)+j]);
            edgesVer[ i       * (nx + 1) + j    ]->neibCells.push_back(cells[i*(nx)+j]);

            cells[i*(nx)+j]->number = i*(nx)+j;
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


void Mesh2D::exportMesh() const
{
    writer << "Lx = " << Lx << endl;
    writer << "Ly = " << Ly << endl;

    writer << "nx = " << nx << endl;
    writer << "ny = " << ny << endl;

    writer << "hx = " << cells[0]->h().x() << endl;
    writer << "hy = " << cells[0]->h().y() << endl;

    writer << "ff = " << nShapes << endl;

    writer << "cellCenters" << endl;

    for (int i = 0; i < nCells; ++i)
    {
        writer << cells[i]->getCellCenter().x() << ' ' << cells[i]->getCellCenter().y() << endl;
    }

    cout << "Mesh export OK" << endl;
}
