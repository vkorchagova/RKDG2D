#include <iostream>
#include "DecomposerUNV.h"
#include <cstdlib>

using namespace std;


// ----------------------------------------------------------------------------


int main()
{

    DecomposerUNV converter("..//..//..//meshCollection//Mesh_forwardStep_005.unv","mesh2D");
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
