#include <iostream>
#include "DecomposerUNV.h"
#include <cstdlib>

using namespace std;


// ----------------------------------------------------------------------------


int main(int argc, char** argv)
{

    DecomposerUNV converter("..//Mesh_1.unv","mesh2D");
    //DecomposerUNV converter("..//..//..//mCol//fs64triag.unv","mesh2D");
    //DecomposerUNV converter("..//..//..//..//tests//unv//dmAngle0p0025triag.unv","mesh2D");
    int nDomains = 0;

    if (argc > 1)
    {
        nDomains = atoi(argv[1]);
    }
    else
    {
        cout << "Number of subdomains is not found." << endl;
        cout << "Please set them as an argument, for example: ./decomposeUNV 2" << endl;
        exit(1);
    }

    if (nDomains < 2)
    {
        cout << "Number of subdomains should be more than 1." << endl;
        exit(1);
    }
    
    string metisCommand = "mpmetis meshMETIS " + to_string(nDomains);
    string partCellsFile = "meshMETIS.epart." + to_string(nDomains);

    converter.importUNV();
    converter.exportRKDG();
    
    converter.exportMETIS();
    system(metisCommand.c_str());

    converter.importPartition(partCellsFile.c_str());
    
    for (int iDom = 0; iDom < nDomains; ++iDom)
        converter.exportPartMeshRKDG(iDom);
    
    converter.exportVTK();
    
    cout << "END" << endl;
    return 0;
}
