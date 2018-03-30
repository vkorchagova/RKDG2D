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

Mesh2D::Mesh2D(string fileName, const Problem& prb)
{
    importMesh(fileName,prb);
}

Mesh2D::~Mesh2D()
{

}

// ------------------ Private class methods --------------------

void Mesh2D::createRectangularMesh(const Problem &prb)
{
    // define number of nodes, edges, cells

//    int nNodes = (nx + 1)*(ny + 1);
//    int nEdgesHor = nx*(ny + 1);
//    int nEdgesVer = ny*(nx + 1);

//    nCells = nx*ny;

//    // space step in x, y direction

//    double hx = Lx / nx;
//    double hy = Ly / ny;

//    // reserve memory

//    nodes.reserve(nNodes);

//    vector<shared_ptr<Edge>> edgesHor;
//    vector<shared_ptr<Edge>> edgesVer;

//    edgesHor.reserve(nEdgesHor);
//    edgesVer.reserve(nEdgesVer);

//    //edgesInternal.reserve( (nx - 1) * (ny - 1) );
//    //edgesBoundary.reserve( 2*nx + 2*ny );
//    edges.reserve( (nx - 1) * (ny - 1) + 2*nx + 2*ny );
//    // fill nodes

//    for (int i = 0; i < ny + 1; ++i)
//        for (int j = 0; j < nx + 1; ++j)
//            nodes.emplace_back( Point({ j * hx, i * hy }) );

//    for (int i = 0; i < nNodes; ++i)
//        nodes[i].number = i;

//    // get horizontal edges

//    for (int j = 0; j < nx; ++j)
//    { // bottom boundary
//        edgesHor.emplace_back( make_shared<EdgeBoundary>(nodes[j], nodes[j + 1]) );
//        edgesHor.back()->n = Point({ 0.0, -1.0 });
//        edgesBoundary.emplace_back(dynamic_pointer_cast<EdgeBoundary>(edgesHor.back()));
//    }


//    for (int i = 1; i < ny; ++i)
//       for (int j = 0; j < nx; ++j)
//       {
//            edgesHor.emplace_back( make_shared<EdgeInternal>(nodes[i*(nx + 1) + j], nodes[i*(nx + 1) + j + 1]) );
//            edgesHor.back()->n = Point({ 0.0, 1.0 });
//            edgesInternal.emplace_back( dynamic_pointer_cast<EdgeInternal>(edgesHor.back()) );
//       }

//    for (int j = 0; j < nx; ++j)
//    { // top boundary
//        edgesHor.emplace_back( make_shared<EdgeBoundary>(nodes[ny*(nx + 1) + j], nodes[ny*(nx + 1) + j + 1]) );
//        edgesHor.back()->n = Point({ 0.0, 1.0 });
//        edgesBoundary.emplace_back( dynamic_pointer_cast<EdgeBoundary>(edgesHor.back()) );
//    }

//    // get vertical edges: boundary + internal

//    for (int i = 0; i < ny; ++i)
//    {
//        // left boundary
//        edgesVer.emplace_back( make_shared<EdgeBoundary>(nodes[i*(nx + 1)], nodes[i*(nx + 1) + nx + 1]) );
//        edgesVer.back()->n = Point({ -1.0, 0.0 });
//        edgesBoundary.emplace_back( dynamic_pointer_cast<EdgeBoundary>(edgesVer.back()) );

//        for (int j = 1; j < nx; ++j)
//        {
//            edgesVer.emplace_back( make_shared<EdgeInternal>(nodes[i*(nx + 1) + j], nodes[i*(nx + 1) + j + nx + 1]) );
//            edgesVer.back()->n = Point({ 1.0, 0.0 });
//            edgesInternal.emplace_back( dynamic_pointer_cast<EdgeInternal>(edgesVer.back()) );
//        }
//        // right boundary
//        edgesVer.emplace_back( make_shared<EdgeBoundary>(nodes[i*(nx + 1) + nx ], nodes[i*(nx + 1) + nx + nx + 1]) );
//        edgesVer.back()->n = Point({ 1.0, 0.0 });
//        //edgesBoundary.emplace_back( dynamic_pointer_cast<EdgeBoundary>(edgesVer.back()) );
//    }

//    // add right boundary to list of boundary cells
//    for (int i = 0; i < nEdgesVer-nx; ++i)
//        if (i % (nx+1) == 0)
//            edgesBoundary.emplace_back( dynamic_pointer_cast<EdgeBoundary>(edgesVer[i+nx]) );

//    //add numbers of edges
//    for (int i = 0; i < edgesBoundary.size(); ++i)
//        edgesBoundary[i]->number = i;

//    for (int i = 0; i < edgesInternal.size(); ++i)
//        edgesInternal[i]->number = i + edgesBoundary.size();


//    // get cells as edges (counter-clockwise: lb -> rb -> ru -> rl)
//    for (int i = 0; i < ny; ++i)
//    {
//        for (int j = 0; j < nx; ++j)
//        {
//            //define list of edges

//            vector<shared_ptr<Edge>> edges = {edgesHor[i*(nx)+j], \
//                                                   edgesVer[(i)*(nx + 1) + j + 1], \
//                                                   edgesHor[(i + 1)*(nx)+j], \
//                                                   edgesVer[(i)*(nx + 1) + j] };

//            // add cells in list
//            cells.emplace_back( make_shared<Cell>(edges,prb));

//            // this cell is neighbour for its edges

//            edgesHor[ i       *  nx      + j    ]->neibCells.push_back(cells[i*(nx)+j]);
//            edgesVer[ i       * (nx + 1) + j + 1]->neibCells.push_back(cells[i*(nx)+j]);
//            edgesHor[(i + 1)  *  nx      + j    ]->neibCells.push_back(cells[i*(nx)+j]);
//            edgesVer[ i       * (nx + 1) + j    ]->neibCells.push_back(cells[i*(nx)+j]);

//            cells[i*(nx)+j]->number = i*(nx)+j;
//        }
//    }

//    // make groups of boundaries
//    patches.resize(4);

//    patches[0].patchName = "bottom";
//    patches[1].patchName = "top";
//    patches[2].patchName = "left";
//    patches[3].patchName = "right";

//    for (int i = 0; i < nx; ++i)
//    {
//        patches[0].edgeGroup.push_back(edgesBoundary[i]);
//        patches[1].edgeGroup.push_back(edgesBoundary[i + nx]);
//    }

//    for (int i = 0; i < ny; ++i)
//    {
//        patches[2].edgeGroup.push_back(edgesBoundary[i + 2*nx]);
//        patches[3].edgeGroup.push_back(edgesBoundary[i + 2*nx + ny]);
//    }
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
        writer << nodes[i].x() << ' ' << nodes[i].y() << endl;

    writer << "$EndNodes\n";

    // export edges             ---------------------------

    writer << "$Edges\n";

    writer << nBoundEdges << endl;
    writer << nEdges << endl;

    for (size_t i = 0; i < nEdges; ++i)
        writer << edges[i]->nodes[0]->number + 1 << ' ' << edges[i]->nodes[1]->number + 1 << endl;

//    for (size_t i = 0; i < edgesInternal.size(); ++i)
//        writer << edgesInternal[i]->nodes[0]->number + 1  << ' ' << edgesInternal[i]->nodes[1]->number + 1 << endl;

    writer << "$EndEdges\n";

    // export cells             ---------------------------

    writer << "$Cells\n";

    writer << cells.size() << endl;

    for (size_t i = 0; i < cells.size(); ++i)
    {
        writer << cells[i]->edges.size();

        for (int j = 0; j < cells[i]->edges.size(); j++)
            writer << ' ' << cells[i]->edges[j]->number  + 1 ;

        writer << endl;
    }

    writer << "$EndCells\n";

    // export cell centers      ---------------------------

    writer << "$CellCenters\n";

    writer << cells.size() << endl;

    for (size_t i = 0; i < cells.size(); ++i)
        writer << cells[i]->getCellCenter().x() << ' ' << cells[i]->getCellCenter().y() << endl;

    writer << "$EndCellCenters\n";

    // export adjoint cells for edges  ---------------------------

    writer << "$AdjointCellsForEdges\n";

    writer << nEdges << endl;

    for (size_t i = 0; i < nEdges; ++i)
    {
        for (size_t j = 0; j < edges[i]->neibCells.size(); ++i)
            writer << edges[i]->neibCells[j]->number + 1 << ' ';
        writer << endl;
    }

    writer << "$EndAdjointCellsForEdges\n";

    // export normals to edges  ---------------------------

    writer << "$EdgeNormals\n";

    writer << nEdges << endl;

    for (size_t i = 0; i < nEdges; ++i)
        writer << edges[i]->n.x() << ' ' << edges[i]->n.y() << endl;

    writer << "$EndEdgeNormals\n";



    // export patches           ---------------------------

    writer << "$Patches\n";

    writer << patches.size() << endl;

    for (size_t i = 0; i < patches.size(); ++i)
    {
        writer << patches[i].patchName << endl;

        writer << patches[i].edgeGroup.size() << endl;

        for (size_t j = 0; j < patches[i].edgeGroup.size(); ++j)
            writer << patches[i].edgeGroup[j]->number + 1  << endl;
    }

    writer << "$EndPatches\n";

    //close writer

    writer.close();
}

void Mesh2D::importMesh(string fileName, const Problem& prb)
{
    string tag;

    reader.open(fileName.c_str());

    if (!reader.is_open())
    {
        cout << "File " << fileName << " is not found\n";
        exit(0);
    };

    cout << "Processing mesh file " << fileName << endl;

    while (reader.peek() != EOF)
    {
        getline(reader, tag);

        if (tag == "$Nodes")
        {
            int nNodes;
            double x, y;

            reader >> nNodes;
            nodes.reserve(nNodes);

            for (int i = 0; i < nNodes; ++i)
            {
                reader >> x >> y;
                nodes.emplace_back(Point({x,y}));
            }

            do
            {
                getline(reader, tag);

            } while (tag != "$EndNodes");

//            for (int i = 0; i < nNodes; ++i)
//                cout << nodes[i].x() << endl;

            cout << "Number of nodes: " << nNodes << endl;
        }
        else if (tag == "$Edges")
        {
            int nBoundEdges;
            int nEdges;
            int node1, node2;

            reader >> nBoundEdges >> nEdges;

            int nInternalEdges = nEdges - nBoundEdges;

            edges.reserve(nEdges);
            //edgesInternal.reserve(nInternalEdges);

            for (int i = 0; i < nBoundEdges; ++i)
            {
                reader >> node1 >> node2;
                edges.emplace_back(make_shared<EdgeBoundary>(nodes[node1-1], nodes[node2-1]));
            }

            for (int i = 0; i < nInternalEdges; ++i)
            {
                reader >> node1 >> node2;
                edges.emplace_back(make_shared<EdgeInternal>(nodes[node1-1], nodes[node2-1]));
            }

//            for (int i = 0; i < nEdges; ++i)
//                cout << edges[i]->nodes[0]->x() << endl;

            do
            {
                getline(reader, tag);

            } while (tag != "$EndEdges");

            cout << "Number of edges: " << nEdges << ", on boundary: " << nBoundEdges << endl;
        }
        else if (tag == "$Cells")
        {
            reader >> nCells;
            cells.reserve(nCells);

            for (int i = 0; i < nCells; ++i)
            {
                int nEdgesInCell;
                int edge;
                reader >> nEdgesInCell;

                cout << nEdgesInCell << endl;

                vector<shared_ptr<Edge>> edges;
                edges.reserve(nEdgesInCell);

                cout << "ur\n";
                for (int j = 0; j < nEdgesInCell; ++j)
                {
                    reader >> edge;
                    cout<< edge << endl;
                    cout << edges[edge-1]->nodes[0]->x() << endl;
                    edges.push_back(edges[edge-1]);
                }

                // add cells in list
                cells.emplace_back( make_shared<Cell>(edges,prb));

                do
                {
                    getline(reader, tag);

                } while (tag != "$EndCells");

            }

            cout << "Number of cells: " << nCells << endl;
        }
        else
        {
            cout << "Tag \"" << tag << "\" is not supported\n";
            exit(0);
            break;
        }
    }

    reader.close();
}
