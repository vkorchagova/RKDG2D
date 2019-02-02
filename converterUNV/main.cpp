#include <iostream>
#include "FileConverter.h"
#include <cstdlib>

using namespace std;


// ----------------------------------------------------------------------------


int main()
{

    FileConverter converter("..//..//..//mCol//Mesh_Sode_1.unv","mesh2D");


    converter.importUNV();
    converter.exportRKDG();
    
    cout << "END" << endl;
    return 0;
}
