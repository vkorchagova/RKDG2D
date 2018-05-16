#include <iostream>
#include "FileConverter.h"

using namespace std;


// ----------------------------------------------------------------------------


int main()
{
    FileConverter converter("..//meshCollection//Mesh_Sode_1.unv","..//mesh2D");

    converter.importUNV();
    converter.exportRKDG();

    cout << "END" << endl;
    return 0;
}
