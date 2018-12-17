#include "DecomposerUNV.h"

#include <functional>
#include <math.h>
#include <ctime>

using namespace std;

// TODO: make it more flexible, may be, via templates!!!!
std::vector<int> parseStringInt(std::string str)
{
    std::vector<int> data;
    int num;

    std::istringstream sstreamer(str);

    while (sstreamer >> num)
        data.push_back(num);

    return data;
}

// TODO: make it more flexible, may be, via templates!!!!
std::vector<double> parseStringDouble(std::string str)
{
    std::vector<double> data;
    double num;

    std::istringstream sstreamer(str);

    while (sstreamer >> num)
        data.push_back(num);

    return data;
}

// -------------------- Constructors & Destructor -----------------------

DecomposerUNV::DecomposerUNV(std::string unvMeshFile, std::string rkdgMeshFile)
{
    reader.open(unvMeshFile.c_str());

    if (!reader.is_open())
    {
        cout << "File " << unvMeshFile << " is not found\n";
        exit(0);
    }

    writer.open(rkdgMeshFile.c_str());

    if (!writer.is_open())
    {
        cout << "File " << rkdgMeshFile << " could not be opened\n";
        exit(0);
    }

    //TODO: check if file format is correct
    
    nNodes = 0;
    nEdges = 0;
}

DecomposerUNV::~DecomposerUNV()
{
    writer.close();
    reader.close();
}

// -------------------------- Private methods ----------------------------

int DecomposerUNV::readTag()
{
    string str;

    do
    {
        getline(reader, str);
    } while (str != SEPARATOR);

    int num;

    reader >> num;
    getline(reader, str);

    return num;
}


void DecomposerUNV::skipSection()
{
    string str;

    do
    {
        getline(reader, str);
    } while (str != SEPARATOR);

    cout << "skip\n";
}


void DecomposerUNV::readNodes()
{
    string str;

    vector<int> nodeChar;
    vector<double> nodeCoord;

    int num;

    double x, y, z;

    do
    {
        getline(reader, str);

        nodeChar = parseStringInt(str);

        if (nodeChar[0] == -1) // end of section
            break;

        getline(reader, str);
        //cout << "rn str: " << str << endl;

        nodeCoord = parseStringDouble(str);
        nodes.push_back(nodeCoord);

    } while (str != SEPARATOR);
    
    nNodes = nodes.size();

    cout << "OK\n";
}  // End readNodes


void DecomposerUNV::readElements()
{

    string str = "";

    vector<int> elementProperties;
    vector<int> beamElementProperties;
    vector<int> elementNodeNumbers;

    int elementType = -100500;

    do
    {
        getline(reader, str);
        //cout << str << endl;

        elementProperties = parseStringInt(str);

        if (elementProperties[0] == -1)
            break;

        elementType = elementProperties[1];

        //cout << "elemType = " << elementType << endl;

        switch (elementType)
        {
            case 11:    // rod
            {
                getline(reader, str);
                beamElementProperties = parseStringInt(str);
                getline(reader, str);
                elementNodeNumbers = parseStringInt(str);

                edges.push_back(elementNodeNumbers[0]);
                edges.push_back(elementNodeNumbers[1]);
                nEdges++;

                break;
            }
            case 41:    // plane triangular
            case 44:    // plane quadrilateral
            {
                getline(reader, str);
                elementNodeNumbers = parseStringInt(str);
                cellsAsNodes.push_back(elementNodeNumbers);

                getElementEdges(elementNodeNumbers);
                getCellCenter(elementNodeNumbers);

                break;
            }
            default:
            {
                cout << "Element type " << elementType << " is not supported\n";
                break;
            }
        }

    } while (str != SEPARATOR);

    nCells = cellsAsEdges.size();
    cout << "OK\n";

}  // End readElements


void DecomposerUNV::readPatches()
{
    string str;

    do
    {
        vector<int> edgeGroup;
        vector<int> patchProperties;

        getline(reader, str);
        patchProperties = parseStringInt(str);

        if (patchProperties[0] == -1)
            break;

        int nEdgesInGroup = patchProperties.back();

        getline(reader, str);

        patchNames.push_back(str);

        //getline(reader, str);

        function<int()> getEdgeNumber = [&]()
        {
            int number;

            for (int j = 0; j < 2; ++j)
                reader >> number;

            int edgeNumber = number;

            for (int j = 0; j < 2; ++j)
                reader >> number;

            return edgeNumber;
        };

        for (int i = 0; i < nEdgesInGroup; ++i)
            edgeGroup.push_back(getEdgeNumber());

        patchEdgeGroups.push_back(edgeGroup);

        getline(reader, str);


    } while (str != SEPARATOR);

    cout << "OK\n";
}  // End readPatches


int DecomposerUNV::checkForExistingEdges (int iNode1, int iNode2) const
{
    //int nExistingEdges = edges.size();
    int pos = 0;
    for (int j = 0; j < nEdges; ++j)
    { // if this edge is exists
        pos = 2 * j;
        if ((iNode1 == edges[pos] && iNode2 == edges[pos + 1]) ||\
            (iNode2 == edges[pos] && iNode1 == edges[pos + 1]))
        {
            return j + 1;
            //elementEdges.push_back(j+1);
            //return true;
        }
    }
    
    return -1;
};

void DecomposerUNV::createNewEdge (int iNode1, int iNode2)
{
    //vector<int> newEdge = {iNode1,iNode2};
    edges.push_back(iNode1);
    edges.push_back(iNode2);
    nEdges += 1;
};


void DecomposerUNV::getElementEdges(const std::vector<int>& nodeNumbers)
{
    int n = nodeNumbers.size();

    vector<int> elementEdges;
    elementEdges.reserve(n);
    

    //- end of useful internal functions
    
    int iEdge = -1;

    for (int i = 0; i < n - 1; ++i)
    {
        iEdge = checkForExistingEdges(nodeNumbers[i],nodeNumbers[i + 1]);
        
        if (iEdge != -1)
        {
            elementEdges.push_back(iEdge);
        }
        else
        {
            createNewEdge(nodeNumbers[i],nodeNumbers[i + 1]);
            elementEdges.push_back(nEdges);
        }
    }
    
    iEdge = checkForExistingEdges(nodeNumbers[n - 1],nodeNumbers[0]);
    
    if (iEdge != -1)
        elementEdges.push_back(iEdge);
    else
    {
        createNewEdge(nodeNumbers[n - 1],nodeNumbers[0]);
        elementEdges.push_back(nEdges);
    }
    
    cellsAsEdges.push_back(elementEdges);

    //cout << "nedgs = " << nEdges << endl;
} // End getElementEdges


void DecomposerUNV::getCellCenter(const vector<int> &nodeNumbers)
{
    vector<double> center = {0.0, 0.0};

    for (int iNode : nodeNumbers )
    {
        center[0] += nodes[iNode - 1][0];
        center[1] += nodes[iNode - 1][1];
    }

    center[0] /= double(nodeNumbers.size());
    center[1] /= double(nodeNumbers.size());

    cellCenters.push_back(center);

} // End getCellCenter


void DecomposerUNV::setAdjointCells()
{
    //int nEdges = edges.size();
    adjEdgeCells.resize(nEdges);
    
    for (size_t iCell = 0; iCell < nCells; ++iCell)
        for (int iEdge : cellsAsEdges[iCell])
            adjEdgeCells[iEdge - 1].push_back(iCell + 1);

} // End setAdjointCells


void DecomposerUNV::getEdgeNormals()
{
    edgeNormals.reserve(edges.size());
    
    vector<double> edgeV(2);
    vector<double> node1(2); // | fast access to nodes of considered edge
    vector<double> node2(2); // |
    vector<double> centerFirstAdjCell(2);
    vector<double> centerV(2);


    for (size_t i = 0; i < nEdges; i++)
    {
        //vector<int> edge = edges[i];
        node1 = nodes[ edges[2*i + 1] - 1];
        node2 = nodes[ edges[2*i] - 1];
        
        edgeV[0] = node1[0] - node2[0];
        edgeV[1] = node1[1] - node2[1];
        
        centerFirstAdjCell = cellCenters[adjEdgeCells[i][0] - 1];
        
        centerV[0] = centerFirstAdjCell[0] - node2[0];
        centerV[1] = centerFirstAdjCell[1] - node2[1];
        
        vector<double> n = { edgeV[1], -edgeV[0]};
        
        if (n[0]*centerV[0] + n[1]*centerV[1] > 0.0)
        {
            n[0] *= -1;
            n[1] *= -1;
        }

        double ilen = 1.0/sqrt(n[0]*n[0] + n[1]*n[1]);

        n[0] *= ilen;
        n[1] *= ilen;

        edgeNormals.push_back(n);
    }
} // End getEdgeNormals


// ------------------------------ Public methods ------------------------------


void DecomposerUNV::importUNV ()
{
    int num = -1;

    while (reader.peek() != EOF)
    {
        num = readTag();

        cout << num << endl;

        switch (num)
        {
            case 164:
            {
                cout << "Processing units ... ";
                skipSection();
                break;
            }
            case 2420:
            {
                cout << "Processing coordinate system ... ";
                skipSection();
                break;
            }
            case 2411:
            {
                cout << "Processing nodes ... ";
                readNodes();
//                for (int i = 0; i < nodes.size(); ++i)
//                    cout << nodes[i][0] << ' ' << nodes[i][1] << endl;
                break;
            }
            case 2412:
            {
                cout << "Read elements (boundary edges + cells) ... ";
                clock_t to,te;
                to = clock();
                readElements();
                te = clock();
                
                cout << "OK" << endl;
                cout << "time = " << float(te - to) / CLOCKS_PER_SEC << endl;
                
                cout << "Processing adjoint cells for each edge ... ";
                
                to = clock();
                setAdjointCells();
                te = clock();
                
                cout << "OK" << endl;
                cout << "time = " << float(te - to) / CLOCKS_PER_SEC << endl;
                
                cout << "Processing normals ... ";
                to = clock();
                getEdgeNormals();
                te = clock();
                
                cout << "OK" << endl;
                cout << "time = " << float(te - to) / CLOCKS_PER_SEC << endl;
                ;

                break;
            }
            case 2467:  //patches
            {
                cout << "Processing patches ... ";
                readPatches();

                cout << "patch names\n";
                for (size_t i = 0; i < patchNames.size(); ++i)
                    cout << i+1 << '|' << patchNames[i] << ", " << patchEdgeGroups[i].size() << " edges" << endl;

                break;
            }
            default:
            {
                cout << "Section " << num << "is not supported\n";
                skipSection();
                break;
            }
        };
    }
}

void DecomposerUNV::exportRKDG()
{
    //renumerateEdges();

    writer.precision(16);

    // --------------------------------------

    writer << "$Nodes\n";
    writer << nNodes << endl;

    for (int i = 0; i < nNodes; ++i)
        writer << nodes[i][0] << ' ' << nodes[i][1] << endl;

    writer << "$EndNodes\n";

    // --------------------------------------

    int nEdgesBound = 0;

    for (int i = 0; i < patchNames.size(); ++i)
        nEdgesBound += patchEdgeGroups.size();

    writer << "$Edges\n";

    writer << nEdgesBound << endl;
    writer << nEdges << endl;

    for (size_t i = 0; i < nEdges; ++i)
        writer << abs (int(adjEdgeCells[i].size() - 2)) << ' ' << edges[2*i] << ' ' << edges[2*i + 1] << endl;

    writer << "$EndEdges\n";

    // --------------------------------------

    writer << "$Cells\n";
    writer << nCells << endl;

    for (size_t i = 0; i < nCells; ++i)
    {
        writer << cellsAsEdges[i].size() << ' ';
        for (size_t j = 0; j < cellsAsNodes[i].size(); ++j)
             writer << cellsAsNodes[i][j] << ' ';
        for (size_t j = 0; j < cellsAsEdges[i].size(); ++j)
             writer << cellsAsEdges[i][j] << ' ';
        writer << endl;
    }

    writer << "$EndCells\n";

    // --------------------------------------

//    writer << "$CellCenters\n";
//    writer << cellsAsEdges.size() << endl;

//    for (size_t i = 0; i < cellsAsEdges.size(); ++i)
//        writer << cellCenters[i][0] << ' ' << cellCenters[i][1] << endl;

//    writer << "$EndCellCenters\n";

    // --------------------------------------

    writer << "$AdjointCellsForEdges\n";
    writer << nEdges << endl;

    for (size_t i = 0; i < nEdges; ++i)
    {
        writer << adjEdgeCells[i].size() << ' ';

        for (size_t j = 0; j < adjEdgeCells[i].size(); ++j)
            writer << adjEdgeCells[i][j] << ' ' ;
        writer << endl;
    }

    writer << "$EndAdjointCellsForEdges\n";

    // --------------------------------------

    writer << "$EdgeNormals\n";
    writer << nEdges << endl;

    for (size_t i = 0; i < nEdges; ++i)
        writer << edgeNormals[i][0] << ' ' << edgeNormals[i][1] << endl;

    writer << "$EndEdgeNormals\n";

    // --------------------------------------

    writer << "$Patches\n";
    writer << patchNames.size() << endl;

    for (size_t i = 0; i < patchEdgeGroups.size(); ++i)
    {
        writer  << patchNames[i] << endl;
        writer  << patchEdgeGroups[i].size() << endl;

        for (size_t j = 0; j < patchEdgeGroups[i].size(); ++j)
            writer << patchEdgeGroups[i][j] << endl;
    }

    writer << "$EndPatches\n";

    cout << "Export OK\n";

}

void DecomposerUNV::exportMETIS() const
{
    ofstream writerMETIS("meshMETIS");
    
    writerMETIS << nCells << endl;
    
    for (size_t i = 0; i < nCells; ++i)
    {
        for (size_t j = 0; j < cellsAsNodes[i].size(); ++j)
             writerMETIS << cellsAsNodes[i][j] << ' ';
        writerMETIS << endl;
    }
    
    writerMETIS.close();
}


void DecomposerUNV::importPartition(std::string metisPartFile)
{
    ifstream metisReader(metisPartFile.c_str());
    
    int buffer = 0;
    
    partition.reserve(nCells);
    
    for (int i = 0; i < nCells; ++i)
    {
        metisReader >> buffer;
        partition.emplace_back(buffer);
    }

    metisReader.close();
}

void DecomposerUNV::exportPartMeshRKDG(int nDom) const
{
    cout << "Exporting domain #" << nDom << "... ";
    
    
    
    // step 0: preparation
    string fileName = "mesh2D." + to_string(nDom);
    ofstream writerPart(fileName.c_str());
    
    vector<bool> isWrittenN;
    vector<bool> isWrittenE;
    vector<int> nodesInDom;
    vector<int> edgesInDom;
    vector<vector<int>> patchGroupsInDom;
    
    vector<int> alienCells;
    vector<int> alienCellsDom;
    
    clock_t ts, te;
    
    ts = clock();
    
    // step 1: read numbers of cells in domain
    vector<int> cellsInDom;
    
    for (int i = 0; i < partition.size(); ++i)
    {
        if (partition[i] == nDom)
            cellsInDom.push_back(i);
    }
    
    
    // step 2: prepare info 
    isWrittenN.resize(nNodes);
    fill(isWrittenN.begin(),isWrittenN.end(),false);
    isWrittenE.resize(nEdges);
    fill(isWrittenE.begin(),isWrittenE.end(),false);
    
    for (const int num : cellsInDom)
    {
        // get coordinates of nodes
        for (const int nNode : cellsAsNodes[num])
            if (!isWrittenN[nNode-1])
            {
                nodesInDom.push_back(nNode-1);
                isWrittenN[nNode-1] = true;
            }

        // get edges
        for (const int nEdge : cellsAsEdges[num])
            if (!isWrittenE[nEdge-1])
            {
                edgesInDom.push_back(nEdge-1);
                isWrittenE[nEdge-1] = true;
            }
    }
    
    //step 3: get patches
    patchGroupsInDom.resize(patchEdgeGroups.size());
    vector<int>::iterator founded;
    
    for (int iPatch = 0; iPatch < patchEdgeGroups.size(); ++iPatch)
    {
        for (int iEdge : patchEdgeGroups[iPatch])
        {
            founded = std::find(edgesInDom.begin(),edgesInDom.end(),iEdge - 1);

            if (founded != edgesInDom.end())
                patchGroupsInDom[iPatch].push_back(iEdge - 1);
        }
    }
    
    int nRealEdges = edgesInDom.size();
    int nRealCells = cellsInDom.size();
    
    
    //step 4: get alien cells
    for (int nEdge : edgesInDom)
    {
        for (int cell : adjEdgeCells[nEdge])
        {
            if (partition[cell-1] != nDom && find(alienCells.begin(), alienCells.end(), cell - 1) == alienCells.end())
            {
                alienCells.push_back(cell - 1);
                alienCellsDom.push_back(partition[cell - 1]);
            }
        }
    }
    
    //step 5: add nodes and edges for alien cells to vectors of nodes and edges
    //for (const int num : alienCells)
    for (int i = 0; i < alienCells.size(); ++i)
    {
        int num = alienCells[i];
        // get coordinates of nodes
        for (const int nNode : cellsAsNodes[num])
        {
            if (!isWrittenN[nNode-1])
            {
                nodesInDom.push_back(nNode-1);
                isWrittenN[nNode-1] = true;
            }
        }

        // get edges
        for (const int nEdge : cellsAsEdges[num])
            if (!isWrittenE[nEdge-1])
            {
                edgesInDom.push_back(nEdge-1);
                isWrittenE[nEdge-1] = true;
            }
    }
    
    
    int nAlienCells = alienCells.size();
    
    te = clock();
    
    cout.precision(8);
    cout << "\n=========\n time = " << (float)(te - ts) / CLOCKS_PER_SEC << endl;
    
    //step last: write file
    
    writerPart.precision(16);
    int globalNum = -1;

    // --------------------------------------

    writerPart << "$Nodes\n";
    writerPart << nodesInDom.size() << endl;

    for (int i = 0; i < nodesInDom.size(); ++i)
    {
        globalNum = nodesInDom[i];
        writerPart << globalNum + 1 << ' ' << nodes[globalNum][0] << ' ' << nodes[globalNum][1] << endl;
    }

    writerPart << "$EndNodes\n";

    // --------------------------------------

    writerPart << "$Edges\n";

    writerPart << nRealEdges << endl;
    writerPart << edgesInDom.size() << endl;

    for (size_t i = 0; i < edgesInDom.size(); ++i)
    {
        globalNum = edgesInDom[i];
        writerPart << globalNum + 1 << ' ' << edges[2*globalNum] << ' ' << edges[2*globalNum + 1] << endl;
    }

    writerPart << "$EndEdges\n";

    // --------------------------------------

    writerPart << "$Cells\n";
    writerPart << nRealCells << endl;

    for (size_t i = 0; i < nRealCells; ++i)
    {
        globalNum = cellsInDom[i];
        writerPart << globalNum + 1 << ' ';
        writerPart << cellsAsNodes[globalNum].size() << ' ';
        for (size_t j = 0; j < cellsAsNodes[globalNum].size(); ++j)
             writerPart << cellsAsNodes[globalNum][j] << ' ';
        for (size_t j = 0; j < cellsAsEdges[globalNum].size(); ++j)
             writerPart << cellsAsEdges[globalNum][j] << ' ';
        writerPart << endl;
    }

    writerPart << "$EndCells\n";

    // --------------------------------------

    writerPart << "$NeibProcCells\n";
    writerPart << nAlienCells << endl;

    for (int i = 0; i < nAlienCells; ++i)
    {
        globalNum = alienCells[i];
        writerPart << globalNum + 1 << ' ';
        writerPart << cellsAsNodes[globalNum].size() << ' ';
        for (size_t j = 0; j < cellsAsNodes[globalNum].size(); ++j)
             writerPart << cellsAsNodes[globalNum][j] << ' ';
        for (size_t j = 0; j < cellsAsEdges[globalNum].size(); ++j)
             writerPart << cellsAsEdges[globalNum][j] << ' ';
             
        writerPart << alienCellsDom[i] << ' ';
        
        writerPart << endl;
    }

    writerPart << "$EndNeibProcCells\n";

    // --------------------------------------

    writerPart << "$AdjointCellsForEdges\n";
    writerPart << nRealEdges << endl;

    for (size_t i = 0; i < nRealEdges; ++i)
    {
        globalNum = edgesInDom[i];
        writerPart << adjEdgeCells[globalNum].size() << ' ';

        for (size_t j = 0; j < adjEdgeCells[globalNum].size(); ++j)
            writerPart << adjEdgeCells[globalNum][j] << ' ' ;
        writerPart << endl;
    }

    writerPart << "$EndAdjointCellsForEdges\n";

    // --------------------------------------

    writerPart << "$EdgeNormals\n";
    writerPart << nRealEdges << endl;

    for (size_t i = 0; i < nRealEdges; ++i)
    {
        globalNum = edgesInDom[i];
        writerPart << edgeNormals[globalNum][0] << ' ' << edgeNormals[globalNum][1] << endl;
    }

    writerPart << "$EndEdgeNormals\n";

    // --------------------------------------

    writerPart << "$Patches\n";
    writerPart << patchNames.size() << endl;

    for (size_t i = 0; i < patchEdgeGroups.size(); ++i)
    {
        writerPart  << patchNames[i] << endl;
        writerPart  << patchGroupsInDom[i].size() << endl;

        for (size_t j = 0; j < patchGroupsInDom[i].size(); ++j)
            writerPart << patchGroupsInDom[i][j] + 1 << endl;
    }

    writerPart << "$EndPatches\n";

    cout << "OK\n";
}


void DecomposerUNV::exportVTK() const
{
    //writer.open("mesh2D.vtk");
    
    ofstream writerVTK("mesh2D.vtk");

    int nNodes = nodes.size();
    int nCells = cellsAsNodes.size();

    writerVTK << "# vtk DataFile Version 2.0" << endl;
    writerVTK << "RKDG 2D data" << endl;
    writerVTK << "ASCII" << endl;

    writerVTK << "DATASET POLYDATA" << endl;

    writerVTK << "POINTS " << nNodes << " float" << endl;

    // write coordinates of nodes
    for (int i = 0; i < nNodes; ++i)
        writerVTK << nodes[i][0] << ' ' << nodes[i][1] << ' ' << "0" << endl;

    //get size of polygon list

    int polySize = 0;

    for (int i = 0; i < nCells; ++i)
        polySize += cellsAsNodes[i].size();

    polySize += nCells;

    writerVTK << "POLYGONS " << nCells << ' ' << polySize << endl;


    // write cells using numbers of nodes
    for (const vector<int> cell : cellsAsNodes)
    {
        writerVTK << cell.size() << ' ';

        for (int node : cell)
            writerVTK << node - 1 << ' ';

        writerVTK << endl;
    }

    writerVTK << "CELL_DATA " << nCells << endl;

    writerVTK << "SCALARS split int" << endl;
    writerVTK << "LOOKUP_TABLE default" << endl;

    for (int i = 0; i < partition.size(); ++i)
        writerVTK << partition[i] << endl;
        
    writerVTK << "SCALARS cell_numeration int" << endl;
    writerVTK << "LOOKUP_TABLE default" << endl;

    for (int i = 0; i < nCells; ++i)
        writerVTK << i << endl;


    cout << "Mesh export OK" << endl;

    writerVTK.close();
}


// EOF