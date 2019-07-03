#include <iostream>
#include "DecomposerUNV.h"
#include <cstdlib>

using namespace std;


// ----------------------------------------------------------------------------


int main()
{

    DecomposerUNV converter("..//square1_192.unv","mesh2D");
    //DecomposerUNV converter("..//..//..//mCol//dmAngle001triag.unv","mesh2D");
    //DecomposerUNV converter("..//..//..//..//tests//unv//dmAngle0p0025triag.unv","mesh2D");
    int nDomains = 2;
    
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
