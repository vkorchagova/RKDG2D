#include <iostream>
#include "DecomposerUNV.h"
#include <cstdlib>

using namespace std;


// ----------------------------------------------------------------------------


int main(int argc, char** argv)
{
    if (argc < 3)
    {
        cout << "There should be 2 arguments: file_name (string) and number_of_domains (int)" << endl;
        return -1;
    }

    DecomposerUNV converter(argv[1],"mesh2D");
    //DecomposerUNV converter("..//..//..//mCol//fs64triag.unv","mesh2D");
    //DecomposerUNV converter("..//..//..//..//tests//unv//dmAngle0p0025triag.unv","mesh2D");
    int nDomains = stoi(argv[2]);
    
   // if (nDomains < 2)
   // {
   //     cout << "Number of subdomains should be more than 1." << endl;
  //      exit(1);
   // }
    
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
