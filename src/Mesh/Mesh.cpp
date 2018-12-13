#include "Mesh.h"
#include <iostream>
#include <fstream>

using namespace std;


// ------------------ Constructors & Destructors ----------------


Mesh::Mesh(std::string& fileName)
{
    importMesh(fileName);
}

Mesh::~Mesh()
{

}

// ------------------ Private class methods --------------------

void Mesh::findNeighbourCells(const shared_ptr<Cell>& cell) const
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


void Mesh::importMesh(string& fileName)
{
    string tag;

    ifstream reader;
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

        int globalNum;

        if (tag == "$Nodes")
        {
            int nNodes;
            double x, y;

            reader >> nNodes;
            nodes.reserve(nNodes);

            for (int i = 0; i < nNodes; ++i)
            {
                reader >> globalNum >> x >> y;
                globalNodeNumber.push_back(globalNum);
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
            int node1, node2;

            reader >> nGhostCells >> nEdges;

            edges.reserve(nEdges);

            bool onBoundary = false;

            for (int i = 0; i < nEdges; ++i)
            {
                reader >> globalNum;
                reader >> onBoundary;
                reader >> node1 >> node2;
                globalEdgeNumber.push_back(globalNum);
                edges.emplace_back(
                            Edge (nodes[localNumber(globalNodeNumber, node1-1)],
                                  nodes[localNumber(globalNodeNumber, node2-1)]) );
            }

            do
            {
                getline(reader, tag);

            } while (tag != "$EndEdges");

            cout << "Number of edges: " << nEdges << ", on boundary: " << nGhostCells << endl;
        }
        else if (tag == "$Cells")
        {
            reader >> nRealCells;
            cells.reserve(nRealCells);
            globalCellNumber.reserve(nRealCells);

            for (int i = 0; i < nRealCells; ++i)
            {
                int nEdgesInCell;
                int entity;

                reader >> globalNum;
                reader >> nEdgesInCell;

                globalCellNumber.push_back(globalNum);

                //cout << nEdgesInCell << endl;
                vector<shared_ptr<Point>> curNodes;
                vector<shared_ptr<Edge>> curEdges;

                curNodes.reserve(nEdgesInCell);
                curEdges.reserve(nEdgesInCell);

                for (int j = 0; j < nEdgesInCell; ++j)
                {
                    reader >> entity;
                    curNodes.push_back(
                                make_shared<Point>(
                                    nodes[localNumber(globalNodeNumber, entity - 1)]) );
                }

                for (int j = 0; j < nEdgesInCell; ++j)
                {
                    reader >> entity;
                    curEdges.push_back(
                                make_shared<Edge>(
                                    edges[localNumber(globalEdgeNumber, entity - 1)] ));
                }

                // add cells in list
                cells.emplace_back( Cell(curNodes,curEdges));
            }

            do
            {
                getline(reader, tag);

            } while (tag != "$EndCells");

            cout << "Number of cells: " << nRealCells << endl;
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
                    edges[i].neibCells.push_back(
                                make_shared<Cell>(cells[localNumber(globalCellNumber, adjCell - 1)] ));
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
                edges[i].n = Point({nx, ny});
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

void Mesh::exportMeshVTK(ostream& writer) const
{
    //writer.open("Mesh.vtk");

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
        polySize += cells[i].nEntities;

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

void Mesh::exportMeshVTK_polyvertices(ostream& writer) const
{
    //writer.open("Mesh.vtk");

    writer << "# vtk DataFile Version 2.0" << endl;
    writer << "RKDG 2D data" << endl;
    writer << "ASCII" << endl;

    writer << "DATASET POLYDATA" << endl;

    writer << "POINTS " << nEntitiesTotal << " float" << endl;

    // write coordinates of nodes
    for (int i = 0; i < nCells; ++i)
        for (int j = 0; j < cells[i].nEntities; ++j)
            writer << cells[i].nodes[j]->x() << ' ' << cells[i].nodes[j]->y() << ' ' << "0" << endl;

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
