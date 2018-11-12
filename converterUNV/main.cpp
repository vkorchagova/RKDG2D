#include <iostream>
#include "FileConverter.h"
#include <cstdlib>

using namespace std;


// ----------------------------------------------------------------------------


int main()
{

    FileConverter converter("../sodRectTriag100.unv","mesh2D");


    converter.importUNV();
    //converter.exportRKDG();
    converter.exportMETIS();
    system("mpmetis meshMETIS 12");

    converter.importPartition("meshMETIS.epart.12");
    converter.exportVTK();
    
    cout << "END" << endl;
    return 0;
}
