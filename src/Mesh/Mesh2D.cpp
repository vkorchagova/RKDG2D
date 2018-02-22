#include "Mesh2D.h"
#include <iostream>

using namespace std;

typedef EdgeBoundary edgeBoundaryT;

// ------------------ Constructors & Destructors ----------------


// nx --- number of cells along x axis
// ny --- number of cells along y axis
Mesh2D::Mesh2D(int nx, int ny, double Lx, double Ly, const Problem &prb) : nx(nx), ny(ny), Lx(Lx), Ly(Ly)
{
    createRectangularMesh(prb);

}

Mesh2D::~Mesh2D()
{

}

// ------------------ Private class methods --------------------

void Mesh2D::createRectangularMesh(const Problem &prb)
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

    vector<shared_ptr<Edge>> edgesHor;
    vector<shared_ptr<Edge>> edgesVer;

    edgesHor.reserve(nEdgesHor);
    edgesVer.reserve(nEdgesVer);

    edgesInternal.reserve( (nx - 1) * (ny - 1) );
    edgesBoundary.reserve( 2*nx + 2*ny );

    // fill nodes

    for (int i = 0; i < ny + 1; ++i)
        for (int j = 0; j < nx + 1; ++j)
            nodes.emplace_back( Point({ j * hx, i * hy }) );

    // get horizontal edges

    for (int j = 0; j < nx; ++j)
    { // bottom boundary
        edgesHor.emplace_back( make_shared<EdgeBoundary>(nodes[j], nodes[j + 1]) );
        edgesHor.back()->n = Point({ 0.0, -1.0 });
        edgesBoundary.emplace_back(dynamic_pointer_cast<EdgeBoundary>(edgesHor.back()));
    }


    for (int i = 1; i < ny; ++i)
       for (int j = 0; j < nx; ++j)
       {
            edgesHor.emplace_back( make_shared<EdgeInternal>(nodes[i*(nx + 1) + j], nodes[i*(nx + 1) + j + 1]) );
            edgesHor.back()->n = Point({ 0.0, 1.0 });
            edgesInternal.emplace_back( dynamic_pointer_cast<EdgeInternal>(edgesHor.back()) );
       }

    for (int j = 0; j < nx; ++j)
    { // top boundary
        edgesHor.emplace_back( make_shared<EdgeBoundary>(nodes[ny*(nx + 1) + j], nodes[ny*(nx + 1) + j + 1]) );
        edgesHor.back()->n = Point({ 0.0, 1.0 });
        edgesBoundary.emplace_back( dynamic_pointer_cast<EdgeBoundary>(edgesHor.back()) );
    }

    // get vertical edges: boundary + internal

    for (int i = 0; i < ny; ++i)
    {
        // left boundary
        edgesVer.emplace_back( make_shared<EdgeBoundary>(nodes[i*(nx + 1)], nodes[i*(nx + 1) + nx + 1]) );
        edgesVer.back()->n = Point({ -1.0, 0.0 });
        edgesBoundary.emplace_back( dynamic_pointer_cast<EdgeBoundary>(edgesVer.back()) );

        for (int j = 1; j < nx; ++j)
        {
            edgesVer.emplace_back( make_shared<EdgeInternal>(nodes[i*(nx + 1) + j], nodes[i*(nx + 1) + j + nx + 1]) );
            edgesVer.back()->n = Point({ 1.0, 0.0 });
            edgesInternal.emplace_back( dynamic_pointer_cast<EdgeInternal>(edgesVer.back()) );
        }
        // right boundary
        edgesVer.emplace_back( make_shared<EdgeBoundary>(nodes[i*(nx + 1) + nx ], nodes[i*(nx + 1) + nx + nx + 1]) );
        edgesVer.back()->n = Point({ 1.0, 0.0 });
        //edgesBoundary.emplace_back( dynamic_pointer_cast<EdgeBoundary>(edgesVer.back()) );
    }

    // add right boundary to list of boundary cells
    for (int i = 0; i < nEdgesVer-nx; ++i)
        if (i % (nx+1) == 0)
            edgesBoundary.emplace_back( dynamic_pointer_cast<EdgeBoundary>(edgesVer[i+nx]) );


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
            cells.emplace_back( make_shared<Cell>(edges,prb));

            // this cell is neighbour for its edges

            edgesHor[ i       *  nx      + j    ]->neibCells.push_back(cells[i*(nx)+j]);
            edgesVer[ i       * (nx + 1) + j + 1]->neibCells.push_back(cells[i*(nx)+j]);
            edgesHor[(i + 1)  *  nx      + j    ]->neibCells.push_back(cells[i*(nx)+j]);
            edgesVer[ i       * (nx + 1) + j    ]->neibCells.push_back(cells[i*(nx)+j]);

            cells[i*(nx)+j]->number = i*(nx)+j;
        }
    }

    // make groups of boundaries
    patches.resize(4);

    patches[0].patchName = "bottom";
    patches[1].patchName = "top";
    patches[2].patchName = "left";
    patches[3].patchName = "right";

    for (int i = 0; i < nx; ++i)
    {
        patches[0].edgeGroup.push_back(edgesBoundary[i]);
        patches[1].edgeGroup.push_back(edgesBoundary[i + nx]);
    }

    for (int i = 0; i < ny; ++i)
    {
        patches[2].edgeGroup.push_back(edgesBoundary[i + 2*nx]);
        patches[3].edgeGroup.push_back(edgesBoundary[i + 2*nx + ny]);
    }
}

// ------------------ Public class methods ---------------------


void Mesh2D::exportMesh() const
{
    writer.open("mesh2D");

    cout << "Lx = " << Lx << endl;
    cout << "Ly = " << Ly << endl;

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

    writer.close();
}
