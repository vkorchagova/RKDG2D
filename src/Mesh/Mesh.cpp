#include "Mesh.h"

#include "compService.h"
#include "defs.h"
#include <iostream>
#include <fstream>

using namespace std;


// ------------------ Constructors & Destructors ----------------


Mesh::Mesh(std::string& fileName, Buffers& buf) : buf(buf)
{
    importMesh(fileName);
    
    //cout << "num of cells = " << cells.size() << endl;
    //cout << "num of real cells = " << nRealCells << endl;
    //nEntitiesTotal = 0;

    for (int i = 0; i < cells.size(); ++i)
    {
        //nEntitiesTotal += cells[i]->nEntities;
        cells[i]->setArea();
        cells[i]->setGaussPoints();
        cells[i]->setJacobian();
        cells[i]->setCellCenter();
        cells[i]->number = i;
    }

    //cout << "---" << endl;
    //cout << "num of cells = " << cells.size() << endl;

    for (int i = 0; i < nRealCells; ++i)
    {
        findNeighbourCells(cells[i]);
        findNeighbourCellsVertex(cells[i]);
    }

    ofstream fileVertices("neib_vertex.txt");

    for (int i = 0; i < nRealCells; ++i)
    {
        fileVertices << "cell #" << cells[i]->number << ": ";

        for (int j = 0; j < cells[i]->neibCellsVertex.size(); ++j)
            fileVertices << cells[i]->neibCellsVertex[j]->number << " ";
        fileVertices << endl;
    }

    fileVertices.close();

    structurizeEdges();

    cout << "Mesh import OK" << endl;

    /// MPI OPERATION

    // get total number of cells
    nCellsGlob = 0;
    MPI_Reduce(
        &nRealCells, 
        &nCellsGlob, 
        1, 
        MPI_INT, 
        MPI_SUM, 
        0, 
        MPI_COMM_WORLD
    );

    // first gather for nRealCells for each proc
    buf.nCellsPerProc.resize(numProcsTotal);
    buf.mpiDisplMesh.resize(numProcsTotal);
    MPI_Gather(
        &nRealCells, 
        1, 
        MPI_INT, 
        &(buf.nCellsPerProc[myRank]), 
        1, 
        MPI_INT, 
        0, 
        MPI_COMM_WORLD
    );

    if (myRank == 0)
    {
        buf.mpiDisplMesh[0] = 0;
        
        for (int i = 1; i < numProcsTotal; ++i)
            buf.mpiDisplMesh[i] = buf.mpiDisplMesh[i-1] + buf.nCellsPerProc[i-1]; 

        // cout << "DISPL" << endl;

        // for (int i = 0; i < numProcsTotal; ++i)
        //     cout << buf.mpiDisplMesh[i] << ' ';
        // cout << endl;
    }
    
    
    // if (false)//(myRank == 0)
    // {
    //     cout << "MYrank = " << myRank << "; ";

    //     //for (int len=0; len<numProcsTotal; ++len)\
    //         cout << nCellsPerProc[len] << " ,,,  ";	
    //     cout << endl;
    // }

    // get global numeration of total mesh cells
    buf.globalMap.resize(nCellsGlob);
    MPI_Gatherv(
        &(globalCellNumber[0]),   /*pointer to beginning of sender*/
        nRealCells,                 /*length of send data*/
        MPI_INT,                     
        &(buf.globalMap[0]),        /*pointer to beginning of receiver*/
        &(buf.nCellsPerProc[0]),    /*lengths of arrays of received data*/
        &(buf.mpiDisplMesh[0]),     /*places of received data in receiver*/
        MPI_INT, 
        0,                          /*number of proc who collect data*/
        MPI_COMM_WORLD
    );
    
    // if (myRank == 0)
    // {
    //     cout << "MYrank = " << myRank << "; MAP: ";
    //     for (int len = 0; len < nCellsGlob; ++len)\
    //         cout << buf.globalMap[len] << " , ";
    //     cout << endl;
    // }


    //cout << "rank = " << myRank << "; total mesh size = " << nCellsGlob << endl;   

}

Mesh::~Mesh()
{

}

// ------------------ Private class methods --------------------

void Mesh::structurizeEdges()
{
    for (const Patch& p : patches )
        for (const shared_ptr<Edge> e : p.edgeGroup)
            edges.erase(find(edges.begin(),edges.end(),e));

    nEdgesBound = 0;
    for (const Patch& p : patches )
    {
        edges.insert(edges.begin(),p.edgeGroup.begin(), p.edgeGroup.end());
        nEdgesBound += p.edgeGroup.size();
    }

    for (int i = 0; i < edges.size(); ++i)
        edges[i]->number = i;
}


void Mesh::findNeighbourCells(const std::shared_ptr<Cell>& cell)
{
    for (const shared_ptr<Edge> edge : cell->edges)
    {
        //cout << edge->neibCells[0]->getArea() << endl;
        //for (const shared_ptr<Cell> c : edge->neibCells)
        //    cout << c->center << "; ";
        //cout << endl;
        

        if (edge->neibCells.size() == 2)
        {
            if (edge->neibCells[0] == cell)
                cell->neibCells.push_back(edge->neibCells[1]);
            else
                cell->neibCells.push_back(edge->neibCells[0]);
        }
    }
}


bool Mesh::moveCell(const std::shared_ptr<Cell>& cell, shared_ptr<Cell>& startCell, vector<shared_ptr<Cell>>& neibCellVertexOneEdge) // searching for neib
{
    shared_ptr<Cell> nextCell = startCell;
    for(auto edge : startCell->edges)
    {

        if (edge->neibCells.size() == 2)
        {
            nextCell = (edge->neibCells[0] == startCell) ? edge->neibCells[1] : edge->neibCells[0];
            if ( cell->hasCommonNode(nextCell) && !(cell->inNeibVertList(nextCell, neibCellVertexOneEdge)) && !(cell->inNeibVertList(nextCell, cell->neibCellsVertex)))
            {
                neibCellVertexOneEdge.push_back(nextCell);
                startCell = nextCell;
                return true;
            }
        }
    }

    return false;
}

bool isCCW(const vector<shared_ptr<Cell>>& neibCellVertexOneEdge, const shared_ptr<Cell>& cell)
{
    if (neibCellVertexOneEdge.size() < 2)
        return true;

    Point v1 = neibCellVertexOneEdge[0]->getCellCenter() - cell->getCellCenter();
    Point v2 = neibCellVertexOneEdge[1]->getCellCenter() - cell->getCellCenter();

    if (v1.x() * v2.y() - v1.y() * v2.x() < 0)
        return false;

    return true;
}


void Mesh::filterNeibCellsVertex(const shared_ptr<Cell>& cell)
{
    vector<shared_ptr<Cell>> f = cell->neibCellsVertex;
    Point centerVectorX = f[1]->getCellCenter() - f[0]->getCellCenter();
    Point centerVectorY;

    for (int iCell = 1; iCell < f.size()-1; ++iCell)
    {
        centerVectorY = f[iCell + 1]->getCellCenter() - f[iCell]->getCellCenter();
        double Cross = centerVectorX.x() * centerVectorY.y() - centerVectorX.y() * centerVectorY.x() ;
        if (Cross <= 0)
        {
            f.erase(f.begin() + iCell);
            iCell -= 1;
            centerVectorX = f[iCell + 1]->getCellCenter() - f[iCell]->getCellCenter();
        }
        else
            centerVectorX = centerVectorY;
    }

    centerVectorY = f[0]->getCellCenter() - f[f.size()-1]->getCellCenter();
    double Cross = centerVectorX.x() * centerVectorY.y() - centerVectorX.y() * centerVectorY.x() ;
    if (Cross <= 0)
        f.erase(f.end());

    cell->neibCellsVertex = f;
}

void Mesh::findNeighbourCellsVertex(const std::shared_ptr<Cell>& cell)
{
    shared_ptr<Edge> startEdge;// = cell->edges[0];
    vector<shared_ptr<Cell>> neibCellVertexOneEdge; // temporary array
    Point centerVertexX;

    for (auto e: cell->edges)
    {
        if (e->neibCells.size() == 2)
        {
            startEdge = e;

            shared_ptr<Cell> startCell = (startEdge->neibCells[0]==cell) ? startEdge->neibCells[1] : startEdge->neibCells[0];
            shared_ptr<Cell> nextCell = startCell;

            shared_ptr<Cell> firstCell = startCell;

            neibCellVertexOneEdge = {cell, startCell};

            bool searchStart = true;

            do
            {
                searchStart = moveCell(cell, startCell, neibCellVertexOneEdge);
            } while (searchStart);

            neibCellVertexOneEdge.erase(neibCellVertexOneEdge.begin());

            bool ccw = isCCW(neibCellVertexOneEdge, cell);

            if (!ccw)
                reverse(neibCellVertexOneEdge.begin(),neibCellVertexOneEdge.end());

            for (auto c: neibCellVertexOneEdge)
            {
                bool iswrite = false;
                for (auto ce: cell->neibCellsVertex)
                {
                    if (c == ce )
                    {
                        iswrite = true;
                        break;
                    }
                }
                if (iswrite == false)
                {
                    cell->neibCellsVertex.push_back(c);
                 //break;
                }
            }



    //âîçâðàùàåìñÿ â ïåðâóþ ÿ÷åéêó

            size_t sz =  neibCellVertexOneEdge.size() + 1;  //cell->neibCellsVertex.size();
    
            searchStart = true;
            nextCell = firstCell;
            startCell = firstCell;
            neibCellVertexOneEdge.clear();
            neibCellVertexOneEdge = {cell, startCell};

            do
            {
                searchStart = moveCell(cell, startCell, neibCellVertexOneEdge);
            } while (searchStart);

            neibCellVertexOneEdge.erase(neibCellVertexOneEdge.begin());

    //for (int i= 0; i< neibCellVertexOneEdge.size(); ++i)
    //  cout << neibCellVertexOneEdge[i]->number << ' ';
    //cout << "finish_2" << endl;

    //for (int i = 0; i< cell->neibCellsVertex.size(); ++i)
    //  cout << "! " << cell->neibCellsVertex[i]->number << ' ';
    //cout << "finish_3" << endl;

            ccw = isCCW(neibCellVertexOneEdge, cell);

            if (!ccw)
                reverse(neibCellVertexOneEdge.begin(), neibCellVertexOneEdge.end());

            for (auto c: neibCellVertexOneEdge)
            {
                bool iswrite = false;
                for (auto ce: cell->neibCellsVertex)
                {
                    if (c == ce )
                    {
                        iswrite = true;
                        break;
                    }
                }
                if (iswrite == false)
                {
                 cell->neibCellsVertex.push_back(c);
                 //break;
                }
            }
        }
    }

    //filterNeibCellsVertex(cell);
}



void Mesh::addToProcPatch(const shared_ptr<Cell>& cell, int numProc)
{
    string pName = "procPatch." + to_string(numProc);

    int iPatch = getPatchByName(procPatches, pName);

    if (iPatch == -1)
    {
        iPatch = procPatches.size();
        procPatches.emplace_back(ProcPatch(pName, numProc));
    }

    procPatches[iPatch].cellGroup.push_back(cell);
}

void Mesh::createPhysicalPatch(const vector<shared_ptr<Edge>>& edgeGroup, const string& pName)
{
    patches.emplace_back(Patch(pName));

    Patch& p = patches.back();

    // for (const shared_ptr<Edge>& e : edgeGroup)
    // {
    //     //shared_ptr<Cell> c = makeGhostCell(e);
    //     p.cellGroup.push_back(makeGhostCell(e));
    // }

    p.edgeGroup = edgeGroup;
}

shared_ptr<Cell> Mesh::makeGhostCell(const shared_ptr<Edge>& e)
{
    vector<shared_ptr<Point>> cNodes;
    vector<shared_ptr<Edge>> cEdges;

    shared_ptr<Cell> realCell = e->neibCells[0];

    Point nodeRef = *e->nodes[0];

    // reflect nodes of real cell with respect to edge
    for (const shared_ptr<Point> node : realCell->nodes)
    {
        //if (node->isEqual(nodeRef) || node->isEqual(*(e->nodes[1])))
        //cout << node << ' ' << node.get() << ' ' << e->nodes[0] << endl;
        //cout << node << ' ' << e->nodes[1] << endl;

        if (node == e->nodes[0] || node == e->nodes[1])
        {

            cNodes.push_back(node);
        }
        else
        {
            Point v = rotate(*node - nodeRef, e->n);
            v.x() *= -1;

            nodes.emplace_back(make_shared<Point>(inverseRotate(v, e->n) + nodeRef));
            cNodes.push_back(nodes.back());

        }
    }

    // construct edges of cell
    for (int i = 0; i < realCell->nEntities-1; ++i)
    {
        vector<shared_ptr<Point>> nodesInEdge = {realCell->nodes[i], realCell->nodes[i + 1]};

        shared_ptr<Edge> iE = make_shared<Edge>(Edge(nodesInEdge));

        if (iE->isEqual(*e))
            cEdges.push_back(e);
        else
        {
            edges.push_back(iE);
            cEdges.push_back(edges.back());
        }
    }// 3/4 edges. The last is coming soon...

    // construct last edge
    vector<shared_ptr<Point>> nodesInEdge = {realCell->nodes[0], realCell->nodes[realCell->nEntities-1]};

    shared_ptr<Edge> iE = make_shared<Edge>(Edge(nodesInEdge));

    if (iE->isEqual(*e))
        cEdges.push_back(e);
    else
    {
        edges.push_back(iE);
        cEdges.push_back(edges.back());
    }

    // construct cell
    cells.push_back(make_shared<Cell>(Cell(cNodes,cEdges)));
    e->neibCells.push_back(cells.back());

    return cells.back();
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
            double x, y;

            reader >> nNodes;
            nodes.reserve(nNodes);

            for (int i = 0; i < nNodes; ++i)
            {
                reader >> globalNum >> x >> y;
                globalNodeNumber.push_back(globalNum - 1);
                nodes.emplace_back(make_shared<Point>(Point({x,y})));
                nodes.back()->gNumber = globalNum - 1;
            }

            do
            {
                getline(reader, tag);

            } while (tag != "$EndNodes");

//            for (int i = 0; i < nNodes; ++i)
//                cout << nodes[i]->x() << endl;

            if (myRank == 0)
                cout << "Number of nodes: " << nNodes << endl;
        }
        else if (tag == "$Edges")
        {
            int node1, node2;

            reader >> nRealEdges >> nTotalEdges;

            edges.reserve(nTotalEdges);

            for (int i = 0; i < nTotalEdges; ++i)
            {
                reader >> globalNum;
                reader >> node1 >> node2;
                globalEdgeNumber.push_back(globalNum-1);
                std::vector<shared_ptr<Point>> nodesInEdge;
                nodesInEdge.push_back( nodes[localNumber(globalNodeNumber, node1 - 1)]);
                nodesInEdge.push_back( nodes[localNumber(globalNodeNumber, node2 - 1)]);
                edges.emplace_back( make_shared<Edge>(Edge(nodesInEdge) ));

                edges.back()->number = i;
            }

            //

            do
            {
                getline(reader, tag);

            } while (tag != "$EndEdges");

            if (myRank == 0)
                cout << "Number of edges: " << nRealEdges << endl;
        }
        else if (tag == "$Cells")
        {
            reader >> nRealCells;
            cells.reserve(nRealCells);
            globalCellNumber.reserve(nRealCells);

/*
            for (int i : globalNodeNumber)
                cout << i << ' ';
	    cout << endl;

            for (int i : globalEdgeNumber)
                cout << i << ' ';
            cout << endl;

*/
            for (int i = 0; i < nRealCells; ++i)
            {
                int nEdgesInCell;
                int entity;

                reader >> globalNum;
                reader >> nEdgesInCell;

                globalCellNumber.push_back(globalNum - 1);

                //cout << nEdgesInCell << endl;
                vector<shared_ptr<Point>> curNodes;
                vector<shared_ptr<Edge>> curEdges;

                curNodes.reserve(nEdgesInCell);
                curEdges.reserve(nEdgesInCell);

                for (int j = 0; j < nEdgesInCell; ++j)
                {
                    reader >> entity;
                    //cout << entity - 1 << ' ' << localNumber(globalNodeNumber, entity - 1) << endl;
                    curNodes.push_back(nodes[localNumber(globalNodeNumber, entity - 1)]);
                }

                for (int j = 0; j < nEdgesInCell; ++j)
                {
                    reader >> entity;
                    curEdges.push_back(edges[localNumber(globalEdgeNumber, entity - 1)] );
                }

                // add cells in list
                cells.emplace_back( make_shared<Cell>(Cell(curNodes,curEdges)));
                //cells.back()->number = i;
            }

            do
            {
                getline(reader, tag);

            } while (tag != "$EndCells");

            if (myRank == 0)
                cout << "Number of cells: " << nRealCells << endl;
        }
        // else if (tag == "$NeibProcCells")
        // {
        //     int nProcCells;

        //     reader >> nProcCells;
        //     cells.reserve(nRealCells + nProcCells);
        //     globalCellNumber.reserve(nRealCells + nProcCells);

        //     for (int i = 0; i < nProcCells; ++i)
        //     {

        //         int nEdgesInCell;
        //         int entity;

        //         reader >> globalNum;
        //         reader >> nEdgesInCell;

        //         globalCellNumber.push_back(globalNum - 1);

        //         //cout << nEdgesInCell << endl;
        //         vector<shared_ptr<Point>> curNodes;
        //         vector<shared_ptr<Edge>> curEdges;

        //         curNodes.reserve(nEdgesInCell);
        //         curEdges.reserve(nEdgesInCell);

        //         for (int j = 0; j < nEdgesInCell; ++j)
        //         {
        //             reader >> entity;
        //             //cout << entity - 1 << ' ' << localNumber(globalNodeNumber, entity - 1) << endl;
        //             curNodes.push_back(nodes[localNumber(globalNodeNumber, entity - 1)] );
        //         }

        //         for (int j = 0; j < nEdgesInCell; ++j)
        //         {
        //             reader >> entity;
        //             curEdges.push_back(edges[localNumber(globalEdgeNumber, entity - 1)] );
        //         }

        //         // add cells in list
        //         cells.push_back(make_shared<Cell>(Cell(curNodes,curEdges)));

        //         // add proc cell to patch
        //         int numProc;
        //         reader >> numProc;

        //         addToProcPatch(cells.back(), numProc);
        //     }

        //     do
        //     {
        //         getline(reader, tag);

        //     } while (tag != "$EndNeibProcCells");

        //     nNeibProcs = procPatches.size();

        //     //cout << "Number of neib procs: " << nNeibProcs << endl;
        // }
        else if (tag == "$NeibProcPatches")
        {
            int totalNumProcCells = 0;
            reader >> nNeibProcs;

            for (int iPatch = 0; iPatch < nNeibProcs; ++iPatch)
            {
                int numProc;
                reader >> numProc;

                int nProcCells;
                reader >> nProcCells;
                totalNumProcCells += nProcCells;

                string pName = "procPatch." + to_string(numProc);
                procPatches.emplace_back(ProcPatch(pName, numProc));

                cells.reserve(nRealCells + totalNumProcCells);
                globalCellNumber.reserve(nRealCells + totalNumProcCells);

                for (int i = 0; i < nProcCells; ++i)
                {
                    int nEdgesInCell;
                    int entity;

                    reader >> globalNum;
                    reader >> nEdgesInCell;

                    globalCellNumber.push_back(globalNum - 1);

                    //cout << nEdgesInCell << endl;
                    vector<shared_ptr<Point>> curNodes;
                    vector<shared_ptr<Edge>> curEdges;

                    curNodes.reserve(nEdgesInCell);
                    curEdges.reserve(nEdgesInCell);

                    for (int j = 0; j < nEdgesInCell; ++j)
                    {
                        reader >> entity;
                        //cout << entity - 1 << ' ' << localNumber(globalNodeNumber, entity - 1) << endl;
                        curNodes.push_back(nodes[localNumber(globalNodeNumber, entity - 1)] );
                    }

                    for (int j = 0; j < nEdgesInCell; ++j)
                    {
                        reader >> entity;
                        curEdges.push_back(edges[localNumber(globalEdgeNumber, entity - 1)] );
                    }

                    // add cells in list
                    cells.push_back(make_shared<Cell>(Cell(curNodes,curEdges)));
                    //cells.back()->number = cells.size() - 1;

                    procPatches[iPatch].cellGroup.push_back(cells.back());
                } // for alien cells

                int nInnerProcCells;
                int iCellGlobal;
                reader >> nInnerProcCells;

                for (int i = 0; i < nInnerProcCells; ++i)
                {
                    reader >> iCellGlobal;
                    procPatches[iPatch].innerCellGroup.push_back(cells[localNumber(globalCellNumber, iCellGlobal - 1)]);
                } // for inner cells
            } // for iPatch

            do
            {
                getline(reader, tag);

            } while (tag != "$EndNeibProcPatches");

            if (myRank == 0)
                cout << "Number of neib procs: " << nNeibProcs << endl;
        }
        else if (tag == "$AdjointCellsForEdges")
        {
            int nAdj;
            int adjCell;
            reader >> nAdj;

            if (nAdj != nRealEdges)
            {
                cout << "Broken file: nAdj = " << nAdj << ", nEdges = " << nRealEdges << endl;
                exit(0);
            }

            for (int i = 0; i < nRealEdges; ++i)
            {
               reader >> nAdj;

               for (int j = 0; j < nAdj; ++j)
               {
                    reader >> adjCell;
                    //cout << adjCell - 1 << ' ' << localNumber(globalCellNumber, adjCell - 1) << endl;
                    edges[i]->neibCells.push_back(cells[localNumber(globalCellNumber, adjCell - 1)] );
               }
            }

            do
            {
                getline(reader, tag);

            } while (tag != "$EndAdjointCellsForEdges");

            if (myRank == 0)
                cout << "Adjoint cells for each edge are obtained" << endl;

        }
        else if (tag == "$EdgeNormals")
        {
            double nx, ny;

            int nAdj;
            reader >> nAdj;

            if (nAdj != nRealEdges)
            {
                cout << "Broken file: nAdj = " << nAdj << ", nEdges = " << nRealEdges << endl;
                exit(0);
            }

            for (int i = 0; i < nRealEdges; ++i)
            {
                reader >> nx >> ny;
                edges[i]->n = Point({nx, ny});
            }

            do
            {
                getline(reader, tag);

            } while (tag != "$EndEdgeNormals");

            if (myRank == 0)
                cout << "Normal vectors for each edge are obtained" << endl;
        }
        else if (tag == "$Patches")
        {
            int nPatches;
            string name;

            reader >> nPatches;

            patches.reserve(nNeibProcs + nPatches);

            for (int iP = 0; iP < nPatches; ++iP)
            {
                reader >> name;

                int nEdgesInPatch;
                reader >> nEdgesInPatch;

                vector<shared_ptr<Edge>> edgeGroup;
                edgeGroup.reserve(nEdgesInPatch);

                int iEdge;

                for (int i = 0; i < nEdgesInPatch; ++i)
                {
                    reader >> iEdge;
                    //cout << iEdge << ' ';
                    //cout << localNumber(globalEdgeNumber, iEdge - 1) << endl;
                    edgeGroup.push_back(edges[localNumber(globalEdgeNumber, iEdge - 1)]);
                }

                createPhysicalPatch(edgeGroup, name);
            }
            //for(int i=0;i<edges.size();++i)\
                cout << edges[i].neibCells.size() << endl;

            if (myRank == 0)
                for (const Patch& p : patches)
                {
                    cout << "Patch name: " << p.name << "; number of edges = " << p.edgeGroup.size() << endl;

                    //for (const shared_ptr<Cell> c : p.cellGroup)
                    //    for (const shared_ptr<Point> n : c->nodes)
                    //        cout << n->x() << ' ' << n->y() << endl;
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



void Mesh::exportMeshVTK_polyvertices(ostream& writer) const
{
    //writer.open("Mesh.vtk");

    writer << "# vtk DataFile Version 2.0" << endl;
    writer << "RKDG 2D data" << endl;
    writer << "ASCII" << endl;

    writer << "DATASET POLYDATA" << endl;

    writer << "POINTS " << nNodes << " float" << endl;

    // write coordinates of nodes
    for (int i = 0; i < nRealCells; ++i)
        for (int j = 0; j < cells[i]->nEntities; ++j)
            writer << cells[i]->nodes[j]->x() << ' ' << cells[i]->nodes[j]->y() << ' ' << "0" << endl;

    //get size of polygon list

    int polySize = 0;

    for (int i = 0; i < nRealCells; ++i)
        polySize += cells[i]->nEntities;

    polySize += nRealCells;

    writer << "POLYGONS " << nRealCells << ' ' << nRealCells +  polySize << endl;


    // write cells using numbers of nodes
    int numPolyVertex = 0;
    for (int i = 0; i < nRealCells; ++i)
    {
        writer << cells[i]->nEntities << ' ';

        for (const shared_ptr<Point> node : cells[i]->nodes)
        {
            writer << numPolyVertex << ' ';
            numPolyVertex++;
        }

        writer << endl;
    }

    //cout << "num(-1) = " << cells[cells.size()-1]->edges[2]->nodes[0]->number << endl;
    //cout << "x(-1) = " << cells[cells.size()-1]->edges[2]->nodes[0]->x() << endl;
    //cout << "&num(-1) = " << &(cells[cells.size()-1]->edges[2]->nodes[0]->number) << endl;


    //writer << "CELL_DATA " << nRealCells << endl;
    //writer << "POINT_DATA " << nNodes << endl;


    //cout << "Mesh export OK" << endl;

    //writer.close();
}
