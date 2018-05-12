#include "FileConverter.h"

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

FileConverter::FileConverter(std::string unvMeshFile, std::string rkdgMeshFile)
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

FileConverter::~FileConverter()
{
    writer.close();
    reader.close();
}

// -------------------------- Private methods ----------------------------

int FileConverter::readTag()
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

void FileConverter::skipSection()
{
    string str;

    do
    {
        getline(reader, str);
    } while (str != SEPARATOR);

    cout << "skip\n";
}



void FileConverter::readNodes()
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

void FileConverter::readElements()
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

        elementType = elementProperties[1];

        if (elementProperties[0] == -1)
            break;

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

void FileConverter::readPatches()
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

void FileConverter::getElementEdges(const std::vector<int>& nodeNumbers)
{
    int n = nodeNumbers.size();

    vector<int> elementEdges;
    elementEdges.reserve(n);

    //- useful internal functions

    function<bool(int,int)> checkForExistingEdges = [&](int iNode1, int iNode2)
    {
        for (size_t j = 0; j < edges.size(); ++j)
        { // if this edge is exists and internal
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


void FileConverter::getCellCenter(const vector<int> &nodeNumbers)
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

void FileConverter::setAdjointCells()
{
    int nEdges = edges.size();

    adjEdgeCells.resize(nEdges);

    for (size_t iCell = 0; iCell < cellsAsEdges.size(); ++iCell)
        for (int iEdge : cellsAsEdges[iCell])
            adjEdgeCells[iEdge-1].push_back(iCell+1);

} // End setAdjointCells

void FileConverter::getEgdeNormals()
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


void FileConverter::importUNV ()
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

//                cout << "edges bound \n";
//                for (size_t i = 0; i < edgesBoundary.size(); ++i)
//                    cout << i+1 << '|' << edgesBoundary[i][0] << ' ' << edgesBoundary[i][1] << endl;

//                cout << "edges internal \n";
//                for (size_t i = 0; i < edgesInternal.size(); ++i)
//                    cout << i+1+edgesBoundary.size() << '|' << edgesInternal[i][0] << ' ' << edgesInternal[i][1] << endl;

//                cout << "cells\n";
//                for (size_t i = 0; i < cellsAsEdges.size(); ++i)
//                {
//                    cout << i+1 << '|';
//                    for (size_t j = 0; j < cellsAsEdges[i].size(); ++j)
//                         cout << cellsAsEdges[i][j] << ' ';
//                    cout << endl;
//                }

//                cout << "cell centers\n";
//                for (size_t i = 0; i < cellsAsEdges.size(); ++i)
//                    cout << i+1 << '|' << cellCenters[i][0] << ' ' << cellCenters[i][1] << endl;

//                cout << "adjoint cells for edges\n";
//                for (size_t i = 0; i < adjEdgeCells.size(); ++i)
//                {
//                    cout << i+1 << '|';
//                    for (size_t j = 0; j < adjEdgeCells[i].size(); ++j)
//                        cout << adjEdgeCells[i][j] << ' ' ;
//                    cout << endl;
//                }

//                cout << "edge normals\n";
//                for (size_t i = 0; i < edgeNormals.size(); ++i)
//                    cout << i+1 << '|' << edgeNormals[i][0] << ' ' << edgeNormals[i][1] << endl;

                break;
            }
            case 2467:  //patches
            {
                cout << "Processing patches ... ";
                readPatches();

                cout << "patch names\n";
                for (size_t i = 0; i < patchNames.size(); ++i)
                    cout << i+1 << '|' << patchNames[i] << ", " << patchEdgeGroups[i].size() << " edges" << endl;

//                cout << "patch groups\n";
//                for (size_t i = 0; i < patchEdgeGroups.size(); ++i)
//                {
//                    cout << i+1 << '|';
//                    for (size_t j = 0; j < patchEdgeGroups[i].size(); ++j)
//                        cout << patchEdgeGroups[i][j] << ' ' ;
//                    cout << endl;
//                }

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

void FileConverter::exportRKDG()
{
    writer.precision(16);

    // --------------------------------------

    writer << "$Nodes\n";
    writer << nodes.size() << endl;

    for (int i = 0; i < nodes.size(); ++i)
        writer << nodes[i][0] << ' ' << nodes[i][1] << endl;

    writer << "$EndNodes\n";

    // --------------------------------------
    int nBoundEdges = 0;

    for (int i = 0; i < patchNames.size(); ++i)
        nBoundEdges += patchEdgeGroups[i].size();

    writer << "$Edges\n";
    writer << nBoundEdges << endl;
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
        writer << adjEdgeCells[i].size() << ' ' ;

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

