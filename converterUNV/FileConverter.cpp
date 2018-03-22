#include "FileConverter.h"

using namespace std;

// -------------------------- Private methods ----------------------------

FileConverter::FileConverter(std::string unvMeshFile, std::string rkdgMeshFile)
{
    reader.open(unvMeshFile.c_str());
    writer.open(rkdgMeshFile.c_str());

    //TODO: check if file format correct
}

FileConverter::~FileConverter()
{
    writer.close();
    reader.close();
}

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
}

void FileConverter::readElements()
{

    string str = "dsfsfg";

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

                break;
            }
            case 44:    // plane quadrilateral
            {
                getline(reader, str);
                elementNodeNumbers = parseStringInt(str);

//                for (int i = 0; i < elementNodeNumbers.size(); ++i)
//                    cout << elementNodeNumbers[i] << ' ' ;
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
