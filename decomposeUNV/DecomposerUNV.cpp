#include "DecomposerUNV.h"

#include <functional>
#include <math.h>

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

    //TODO: check if file format correct
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

                edges.push_back(elementNodeNumbers);

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


void DecomposerUNV::getElementEdges(const std::vector<int>& nodeNumbers)
{
    int n = nodeNumbers.size();

    vector<int> elementEdges;
    elementEdges.reserve(n);

    //- useful internal functions

    function<bool(int,int)> checkForExistingEdges = [&](int iNode1, int iNode2)
    {

        int nExistingEdges = edges.size();

        for (int j = 0; j < nExistingEdges; ++j)
        { // if this edge is exists
            if ((iNode1 == edges[j][0] && iNode2 == edges[j][1]) ||\
                (iNode1 == edges[j][1] && iNode2 == edges[j][0]))
            {
                elementEdges.push_back(j+1);
                return true;
            }
        }

        return false;
    };

    function<void(int,int)> createNewEdge = [&](int iNode1, int iNode2)
    {
        vector<int> newEdge = {iNode1,iNode2};
        edges.push_back(newEdge);
        elementEdges.push_back(edges.size());
    };

    //- end of useful internal functions

    for (int i = 0; i < n-1; ++i)
    {
        if (!checkForExistingEdges(nodeNumbers[i],nodeNumbers[i+1]))
        {
           createNewEdge(nodeNumbers[i],nodeNumbers[i+1]);
        }
    }

    if (!checkForExistingEdges(nodeNumbers[n-1],nodeNumbers[0]))
    {
       createNewEdge(nodeNumbers[n-1],nodeNumbers[0]);
    }


    cellsAsEdges.push_back(elementEdges);

} // End getElementEdges


void DecomposerUNV::getCellCenter(const vector<int> &nodeNumbers)
{
    vector<double> center = {0.0, 0.0};

    for (int iNode : nodeNumbers )
    {
        center[0] += nodes[iNode-1][0];
        center[1] += nodes[iNode-1][1];
    }

    center[0] *= 1.0/nodeNumbers.size();
    center[1] *= 1.0/nodeNumbers.size();

    cellCenters.push_back(center);

} // End getCellCenter


void DecomposerUNV::setAdjointCells()
{
    int nEdges = edges.size();

    adjEdgeCells.resize(nEdges);

    for (size_t iCell = 0; iCell < cellsAsEdges.size(); ++iCell)
        for (int iEdge : cellsAsEdges[iCell])
            adjEdgeCells[iEdge-1].push_back(iCell+1);

} // End setAdjointCells


void DecomposerUNV::getEgdeNormals()
{
    edgeNormals.reserve(edges.size());

    for (size_t i = 0; i < edges.size(); ++i)
    {
        vector<int> edge = edges[i];
        vector<double> edgeV = {nodes[edge[1]-1][0] - nodes[edge[0]-1][0], nodes[edge[1]-1][1] - nodes[edge[0]-1][1]};
        vector<double> centerV = {cellCenters[adjEdgeCells[i][0]-1][0] - nodes[edge[0]-1][0], cellCenters[adjEdgeCells[i][0]-1][1] - nodes[edge[0]-1][1]};
        vector<double> n = {edgeV[1], -edgeV[0]};
        if (n[0]*centerV[0] + n[1]*centerV[1] > 0)
        {
            n[0] *= -1;
            n[1] *= -1;
        }

        double ilen = 1.0/sqrt(n[0]*n[0] + n[1]*n[1]);

        n[0] *= ilen;
        n[1] *= ilen;

        edgeNormals.push_back(n);
    }
} // End getEgdeNormals


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
                cout << "Processing elements (boundary edges + cells) ... ";
                readElements();
                setAdjointCells();
                getEgdeNormals();

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
    writer << nodes.size() << endl;

    for (int i = 0; i < nodes.size(); ++i)
        writer << nodes[i][0] << ' ' << nodes[i][1] << endl;

    writer << "$EndNodes\n";

    // --------------------------------------

    int nEdgesBound = 0;

    for (int i = 0; i < patchNames.size(); ++i)
        nEdgesBound += patchEdgeGroups.size();

    writer << "$Edges\n";

    writer << nEdgesBound << endl;
    writer << edges.size() << endl;

    for (size_t i = 0; i < edges.size(); ++i)
        writer << abs (int(adjEdgeCells[i].size() - 2)) << ' ' << edges[i][0] << ' ' << edges[i][1] << endl;

    writer << "$EndEdges\n";

    // --------------------------------------

    writer << "$Cells\n";
    writer << cellsAsEdges.size() << endl;

    for (size_t i = 0; i < cellsAsEdges.size(); ++i)
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
    writer << edges.size() << endl;

    for (size_t i = 0; i < adjEdgeCells.size(); ++i)
    {
        writer << adjEdgeCells[i].size() << ' ';

        for (size_t j = 0; j < adjEdgeCells[i].size(); ++j)
            writer << adjEdgeCells[i][j] << ' ' ;
        writer << endl;
    }

    writer << "$EndAdjointCellsForEdges\n";

    // --------------------------------------

    writer << "$EdgeNormals\n";
    writer << edges.size() << endl;

    for (size_t i = 0; i < edgeNormals.size(); ++i)
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
    int n = cellsAsNodes.size();
    writerMETIS << n << endl;
    for (size_t i = 0; i < n; ++i)
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
    
    partition.reserve(cellsAsNodes.size());
    
    //while (metisReader.peek() != EOF)
    for (int i = 0; i < cellsAsNodes.size(); ++i)
    {
        metisReader >> buffer;
        partition.emplace_back(buffer);
    }

    metisReader.close();
}

void DecomposerUNV::exportPartMeshRKDG(int numDom) const
{
    cout << "Exporting domain #" << numDom << "... ";
    
    // step 0: preparation
    string fileName = "mesh2D." + to_string(numDom);
    ofstream writerPart(fileName.c_str());
    
    vector<bool> isWritten;
    vector<int> nodesInDom;
    vector<int> edgesInDom;
    vector<vector<int>> patchGroupsInDom;
    
    // step 1: read numbers of cells in domain
    vector<int> cellsInDom;
    
    for (int i = 0; i < partition.size(); ++i)
    {
        if (partition[i] == numDom)
            cellsInDom.push_back(i);
    }
    
    // step 2: get coordinates of nodes
    isWritten.resize(nodes.size());
    fill(isWritten.begin(),isWritten.end(),false);
    
    for (const int num : cellsInDom)
        for (const int nNode : cellsAsNodes[num])
            if (!isWritten[nNode])
            {
                nodesInDom.push_back(nNode-1);
                isWritten[nNode] = true;
            }
    
    //step 3: get edges: global_num node_1_global node_2_global
    isWritten.resize(edges.size());
    fill(isWritten.begin(),isWritten.end(),false);
    
    for (const int num : cellsInDom)
        for (const int nEdge : cellsAsEdges[num])
            if (!isWritten[nEdge])
            {
                edgesInDom.push_back(nEdge-1);
                isWritten[nEdge] = true;
            }
            
    //step 4: get patches
    patchGroupsInDom.resize(patchEdgeGroups.size());
    vector<int>::iterator founded;
    
    for (int iPatch = 0; iPatch < patchEdgeGroups.size(); ++iPatch)
    {
        for (int iEdge : patchEdgeGroups[iPatch])
        {
            founded = std::find(edgesInDom.begin(),edgesInDom.end(),iEdge);

            if (founded != edgesInDom.end())
                patchGroupsInDom[iPatch].push_back(iEdge-1);
        }
    }
    
    //step 5: write file
    
    writerPart.precision(16);
    int globalNum = -1;

    // --------------------------------------

    writerPart << "$Nodes\n";
    writerPart << nodesInDom.size() << endl;

    for (int i = 0; i < nodesInDom.size(); ++i)
    {
        globalNum = nodesInDom[i];
        writerPart << globalNum << ' ' << nodes[globalNum][0] << ' ' << nodes[globalNum][1] << endl;
    }

    writerPart << "$EndNodes\n";

    // --------------------------------------

    int nEdgesBound = 0;

    for (int i = 0; i < patchNames.size(); ++i)
        nEdgesBound += patchGroupsInDom.size();

    writerPart << "$Edges\n";

    writerPart << nEdgesBound << endl;
    writerPart << edgesInDom.size() << endl;

    for (size_t i = 0; i < edgesInDom.size(); ++i)
    {
        globalNum = edgesInDom[i];
        writerPart << abs (int(adjEdgeCells[globalNum].size() - 2)) << ' ' << edges[globalNum][0] << ' ' << edges[globalNum][1] << endl;
    }

    writerPart << "$EndEdges\n";

    // --------------------------------------

    writerPart << "$Cells\n";
    writerPart << cellsInDom.size() << endl;

    for (size_t i = 0; i < cellsInDom.size(); ++i)
    {
        globalNum = cellsInDom[i];
        writerPart << globalNum << ' ';
        writerPart << cellsInDom.size() << ' ';
        for (size_t j = 0; j < cellsAsNodes[globalNum].size(); ++j)
             writerPart << cellsAsNodes[globalNum][j] << ' ';
        for (size_t j = 0; j < cellsAsEdges[globalNum].size(); ++j)
             writerPart << cellsAsEdges[globalNum][j] << ' ';
        writerPart << endl;
    }

    writerPart << "$EndCells\n";

    // --------------------------------------

    writerPart << "$AdjointCellsForEdges\n";
    writerPart << edgesInDom.size() << endl;

    for (size_t i = 0; i < edgesInDom.size(); ++i)
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
    writerPart << edgesInDom.size() << endl;

    for (size_t i = 0; i < edgesInDom.size(); ++i)
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
            writerPart << patchGroupsInDom[i][j] << endl;
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


    cout << "Mesh export OK" << endl;

    writerVTK.close();
}


// EOF