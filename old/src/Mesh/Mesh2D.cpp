#include "Mesh2D.h"
#include <iostream>

using namespace std;

typedef EdgeBoundary edgeBoundaryT;

// ------------------ Constructors & Destructors ----------------


Mesh2D::Mesh2D(string fileName, const Problem& prb)
{
    importMesh(fileName,prb);

    nEntitiesTotal = 0;

    for (int i = 0; i < nCells; ++i)
    {
        nEntitiesTotal += cells[i]->nEntities;
        cells[i]->number = i;
        cells[i]->setArea();
        cells[i]->setGaussPoints();
        cells[i]->setJacobian();
        cells[i]->setCellCenter();
        cells[i]->setBasisFunctions();
        cells[i]->setGramian();
        findNeighbourCells(cells[i]);
    }
}

Mesh2D::~Mesh2D()
{

}

// ------------------ Private class methods --------------------

void Mesh2D::findNeighbourCells(const shared_ptr<Cell>& cell) const
{
    for (const shared_ptr<Edge> edge : cell->edges)
    {
        if (edge->neibCells.size() == 2)
        {
            if (edge->neibCells[0] == cell)
                cell->neibCells.push_back(edge->neibCells[1]);
            else
                cell->neibCells.push_back(edge->neibCells[0]);
        }
    }
}

// ------------------ Public class methods ---------------------



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

//        for (int j = 0; j < cells[i]->edges.size(); j++)
//            writer << ' ' << cells[i]->nodes[j]->number  + 1 ;

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

} // end exportMesh

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
                nodes.emplace_back(Node(Point({x,y}),i));
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
            int node1, node2;

            reader >> nBoundEdges >> nEdges;

            edges.reserve(nEdges);
            bool onBoundary = false;

            for (int i = 0; i < nEdges; ++i)
            {
                reader >> onBoundary;
                reader >> node1 >> node2;

                if (onBoundary)
                    edges.emplace_back(make_shared<EdgeBoundary>(nodes[node1-1], nodes[node2-1]));
                else
                    edges.emplace_back(make_shared<EdgeInternal>(nodes[node1-1], nodes[node2-1]));
            }


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
                int entity;
                reader >> nEdgesInCell;

                //cout << nEdgesInCell << endl;
                vector<shared_ptr<Node>> curNodes;
                vector<shared_ptr<Edge>> curEdges;

                curNodes.reserve(nEdgesInCell);
                curEdges.reserve(nEdgesInCell);

                for (int j = 0; j < nEdgesInCell; ++j)
                {
                    reader >> entity;
                    curNodes.push_back(make_shared<Node>(nodes[entity-1]));
                }

                for (int j = 0; j < nEdgesInCell; ++j)
                {
                    reader >> entity;
                    curEdges.push_back(edges[entity-1]);
                }

                // add cells in list
                cells.emplace_back( make_shared<Cell>(curNodes,curEdges,prb));
            }

            do
            {
                getline(reader, tag);

            } while (tag != "$EndCells");

            cout << "Number of cells: " << nCells << endl;
        }
        else if (tag == "$CellCenters")
        {
            int nCenters;
            reader >> nCenters;

            if (nCenters != nCells)
            {
                cout << "Broken file: nCenters != nCells\n";
                exit(0);
            }

            do
            {
                getline(reader, tag);

            } while (tag != "$EndCellCenters");

            //cout << "Cell centers proceeded" << endl;
        }
        else if (tag == "$AdjointCellsForEdges")
        {
            int nAdj;
            int adjCell;
            reader >> nAdj;

            if (nAdj != nEdges)
            {
                cout << "Broken file: nAdj = " << nAdj << ", nEdges = " << nEdges << endl;
                exit(0);
            }

            for (int i = 0; i < nEdges; ++i)
            {
               reader >> nAdj;

               for (int j = 0; j < nAdj; ++j)
               {
                    reader >> adjCell;
                    edges[i]->neibCells.push_back(cells[adjCell - 1]);
               }
            }

            do
            {
                getline(reader, tag);

            } while (tag != "$EndAdjointCellsForEdges");

            //cout << "Adjoint cells for each edge are proceeded" << endl;

        }
        else if (tag == "$EdgeNormals")
        {
            getline(reader, tag);
            double nx, ny;

            for (int i = 0; i < nEdges; ++i)
            {
                reader >> nx >> ny;
                edges[i]->n = Point({nx, ny});
            }

            do
            {
                getline(reader, tag);

            } while (tag != "$EndEdgeNormals");

            //cout << "Normal vectors for each edge are proceeded" << endl;
        }
        else if (tag == "$Patches")
        {
            int nPatches;
            reader >> nPatches;

            patches.resize(nPatches);

            for (int iP = 0; iP < nPatches; ++iP)
            {
                reader >> patches[iP].patchName;

                int nEdgesInPatch;
                reader >> nEdgesInPatch;
                int edge;

                for (int i = 0; i < nEdgesInPatch; ++i)
                {
                    reader >> edge;
                    patches[iP].edgeGroup.push_back(dynamic_pointer_cast<EdgeBoundary>(edges[edge-1]));
                }

                cout << "Patch #" << iP << ": " << patches[iP].patchName << ", " << nEdgesInPatch << " edges contains\n";
            }

            do
            {
                getline(reader, tag);

            } while (tag != "$EndPatches");

           // cout << "Patches are proceeded" << endl;
        }
        else
        {
            cout << "Tag \"" << tag << "\" is not supported\n";
            exit(0);
            break;
        }
    }

    reader.close();
} // end importMesh

void Mesh2D::exportMeshVTK(ostream& writer) const
{
    //writer.open("mesh2D.vtk");

    int nNodes = nodes.size();

    writer << "# vtk DataFile Version 2.0" << endl;
    writer << "RKDG 2D data" << endl;
    writer << "ASCII" << endl;

    writer << "DATASET POLYDATA" << endl;

    writer << "POINTS " << nNodes << " float" << endl;

    // write coordinates of nodes
    for (int i = 0; i < nNodes; ++i)
        writer << nodes[i].x() << ' ' << nodes[i].y() << ' ' << "0" << endl;

    //get size of polygon list

    int polySize = 0;

    for (int i = 0; i < nCells; ++i)
        polySize += cells[i]->nEntities;

    polySize += nCells;

    writer << "POLYGONS " << nCells << ' ' << polySize << endl;


    // write cells using numbers of nodes
    for (const shared_ptr<Cell> cell : cells)
    {
        writer << cell->nEntities << ' ';

        for (const shared_ptr<Node> node : cell->nodes)
            writer << node->number << ' ';

        writer << endl;
    }

    //cout << "num(-1) = " << cells[cells.size()-1]->edges[2]->nodes[0]->number << endl;
    //cout << "x(-1) = " << cells[cells.size()-1]->edges[2]->nodes[0]->x() << endl;
    //cout << "&num(-1) = " << &(cells[cells.size()-1]->edges[2]->nodes[0]->number) << endl;


    //writer << "CELL_DATA " << nCells << endl;
    //writer << "POINT_DATA " << nNodes << endl;


    //cout << "Mesh export OK" << endl;

    //writer.close();
}

void Mesh2D::exportMeshVTK_polyvertices(ostream& writer) const
{
    //writer.open("mesh2D.vtk");

    writer << "# vtk DataFile Version 2.0" << endl;
    writer << "RKDG 2D data" << endl;
    writer << "ASCII" << endl;

    writer << "DATASET POLYDATA" << endl;

    writer << "POINTS " << nEntitiesTotal << " float" << endl;

    // write coordinates of nodes
    for (int i = 0; i < nCells; ++i)
        for (int j = 0; j < cells[i]->nEntities; ++j)
            writer << cells[i]->nodes[j]->x() << ' ' << cells[i]->nodes[j]->y() << ' ' << "0" << endl;

    //get size of polygon list

    writer << "POLYGONS " << nCells << ' ' << nCells +  nEntitiesTotal << endl;


    // write cells using numbers of nodes
    int numPolyVertex = 0;
    for (const shared_ptr<Cell> cell : cells)
    {
        writer << cell->nEntities << ' ';

        for (const shared_ptr<Node> node : cell->nodes)
        {
            writer << numPolyVertex << ' ';
            numPolyVertex++;
        }

        writer << endl;
    }

    //cout << "num(-1) = " << cells[cells.size()-1]->edges[2]->nodes[0]->number << endl;
    //cout << "x(-1) = " << cells[cells.size()-1]->edges[2]->nodes[0]->x() << endl;
    //cout << "&num(-1) = " << &(cells[cells.size()-1]->edges[2]->nodes[0]->number) << endl;


    //writer << "CELL_DATA " << nCells << endl;
    //writer << "POINT_DATA " << nNodes << endl;


    //cout << "Mesh export OK" << endl;

    //writer.close();
}
