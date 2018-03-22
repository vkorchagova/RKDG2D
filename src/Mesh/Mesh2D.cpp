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

Mesh2D::Mesh2D(string fileName)
{
    importMesh(fileName);
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

    for (int i = 0; i < nNodes; ++i)
        nodes[i].number = i;

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

    //add numbers of edges
    for (int i = 0; i < edgesBoundary.size(); ++i)
        edgesBoundary[i]->number = i;

    for (int i = 0; i < edgesInternal.size(); ++i)
        edgesInternal[i]->number = i + edgesBoundary.size();


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


void Mesh2D::exportUniformMesh() const
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


//-----------------------------------------------------------------------

void Mesh2D::exportMesh() const
{
    // open writer

    writer.open("mesh2D");

    // export mesh nodes        ---------------------------

    writer << "$Nodes\n";

    writer << nodes.size() << endl;

    for (size_t i = 0; i < nodes.size(); ++i)
        writer << i+1 << ' ' << nodes[i].x() << ' ' << nodes[i].y() << endl;

    writer << "$EndNodes\n";

    // export edges             ---------------------------

    writer << "$Edges\n";

    writer << edgesBoundary.size() << endl;
    writer << edgesBoundary.size() + edgesInternal.size() << endl;

    for (size_t i = 0; i < edgesBoundary.size(); ++i)
        writer << i+1 << ' ' << edgesBoundary[i]->nodes[0]->number + 1 << ' ' << edgesBoundary[i]->nodes[1]->number + 1 << endl;

    for (size_t i = 0; i < edgesInternal.size(); ++i)
        writer << i+1+edgesBoundary.size() << ' ' << edgesInternal[i]->nodes[0]->number + 1  << ' ' << edgesInternal[i]->nodes[1]->number + 1 << endl;

    writer << "$EndEdges\n";

    // export normals to edges  ---------------------------

    writer << "$EdgeNormals\n";

    writer << edgesBoundary.size() + edgesInternal.size() << endl;

    for (size_t i = 0; i < edgesBoundary.size(); ++i)
        writer << i+1 << ' ' << edgesBoundary[i]->n.x() << ' ' << edgesBoundary[i]->n.y() << endl;

    for (size_t i = 0; i < edgesInternal.size(); ++i)
        writer << i+1+edgesBoundary.size() << ' ' << edgesInternal[i]->n.x() << ' ' << edgesInternal[i]->n.y() << endl;


    writer << "$EndEdgeNormals\n";

    // export cells             ---------------------------

    writer << "$Cells\n";

    writer << cells.size() << endl;

    for (size_t i = 0; i < cells.size(); ++i)
    {
        writer << i+1 << ' ' << cells[i]->edges.size();

        for (int j = 0; j < cells[i]->edges.size(); j++)
            writer << ' ' << cells[i]->edges[j]->number  + 1 ;

        writer << endl;
    }

    writer << "$EndCells\n";

    // export cell centers      ---------------------------

    writer << "$CellCenters\n";

    writer << cells.size() << endl;

    for (size_t i = 0; i < cells.size(); ++i)
        writer << i+1 << ' ' << cells[i]->getCellCenter().x() << ' ' << cells[i]->getCellCenter().y() << endl;

    writer << "$EndCellCenters\n";

    // export patches           ---------------------------

    writer << "$Patches\n";

    writer << patches.size() << endl;

    for (size_t i = 0; i < patches.size(); ++i)
    {
        writer << i + 1  << ' ' << patches[i].patchName << endl;

        writer << patches[i].edgeGroup.size() << endl;

        for (size_t j = 0; j < patches[i].edgeGroup.size(); ++j)
            writer << patches[i].edgeGroup[j]->number + 1  << endl;
    }

    writer << "$EndPatches\n";

    //close writer

    writer.close();
}

void Mesh2D::importMesh(string fileName) const
{
    string x;
    reader.open(fileName);
    getline(reader,x);

    cout << x << endl;

    //first blocks are standard (Cartesian coordinate system, SI units, no additional info about FE)

    //find block with mesh numbers



    reader.close();
}
