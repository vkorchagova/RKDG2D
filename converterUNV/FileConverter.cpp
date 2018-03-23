#include "FileConverter.h"

#include <functional>

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

                edgesBoundary.push_back(elementNodeNumbers);

                break;
            }
            case 44:    // plane quadrilateral
            {
                getline(reader, str);
                elementNodeNumbers = parseStringInt(str);

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

void FileConverter::getElementEdges(const std::vector<int>& nodeNumbers)
{
    int n = nodeNumbers.size();

    vector<int> elementEdges;
    elementEdges.reserve(n);

    //- useful internal functions

    function<bool(int,int)> checkForBoundaryEdges = [&](int iNode1, int iNode2)
    {
        int nBoundEdges = edgesBoundary.size();

        for (int j = 0; j < nBoundEdges; ++j)
        { // if this edge is exists and boundary
            if ((iNode1 == edgesBoundary[j][0] && iNode2 == edgesBoundary[j][1]) ||\
                (iNode1 == edgesBoundary[j][1] && iNode2 == edgesBoundary[j][0]))
            {
                elementEdges.push_back(j+1);
                return true;
            }
        }

        return false;
    };

    function<bool(int,int)> checkForInternalEdges = [&](int iNode1, int iNode2)
    {
        int nInternalEdges = edgesInternal.size();

        for (int j = 0; j < nInternalEdges; ++j)
        { // if this edge is exists and internal
            if ((iNode1 == edgesInternal[j][0] && iNode2 == edgesInternal[j][1]) ||\
                (iNode1 == edgesInternal[j][1] && iNode2 == edgesInternal[j][0]))
            {
                elementEdges.push_back(j+edgesBoundary.size()+1);
                return true;
            }
        }
        return false;
    };

    function<void(int,int)> createNewEdge = [&](int iNode1, int iNode2)
    {
        vector<int> newEdge = {iNode1,iNode2};
        edgesInternal.push_back(newEdge);
        elementEdges.push_back(edgesBoundary.size() + edgesInternal.size());
    };

    //- end of useful internal functions

    for (int i = 0; i < n-1; ++i)
    {
        if (!(checkForBoundaryEdges(nodeNumbers[i],nodeNumbers[i+1]) || \
              checkForInternalEdges(nodeNumbers[i],nodeNumbers[i+1])))
        {
           createNewEdge(nodeNumbers[i],nodeNumbers[i+1]);
        }
    }

    if (!(checkForBoundaryEdges(nodeNumbers[n-1],nodeNumbers[0]) || \
          checkForInternalEdges(nodeNumbers[n-1],nodeNumbers[0])))
    {
       createNewEdge(nodeNumbers[n-1],nodeNumbers[0]);
    }

    cells.push_back(elementEdges);

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

void FileConverter::getEgdeNormals()
{

}


// ------------------------------ Public methods ------------------------------


void FileConverter::importUNV ()
{
    int num = -1;

    while (reader.peek() != EOF)
    {
        num = readTag();

//        cout << num << endl;

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

                cout << "edges bound \n";
                for (size_t i = 0; i < edgesBoundary.size(); ++i)
                    cout << edgesBoundary[i][0] << ' ' << edgesBoundary[i][1] << endl;

                cout << "edges internal \n";
                for (size_t i = 0; i < edgesInternal.size(); ++i)
                    cout << edgesInternal[i][0] << ' ' << edgesInternal[i][1] << endl;

                cout << "cells\n";
                for (size_t i = 0; i < cells.size(); ++i)
                {
                    for (size_t j = 0; j < cells[i].size(); ++j)
                         cout << cells[i][j] << ' ';
                    cout << endl;
                }

                cout << "cell centers\n";
                for (size_t i = 0; i < cells.size(); ++i)
                    cout << cellCenters[i][0] << ' ' << cellCenters[i][1] << endl;

                break;
            }
            case 2467:  //patches
            {
                cout << "Processing patches ... ";
                skipSection();
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


