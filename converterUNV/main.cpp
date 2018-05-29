#include <iostream>
#include "FileConverter.h"

using namespace std;


// ----------------------------------------------------------------------------


int main()
{
    FileConverter converter("..//meshCollection//Mesh_forwardStep_00125.unv","..//mesh2D");

    converter.importUNV();
    converter.exportRKDG();

    cout << "END" << endl;
    return 0;
}
