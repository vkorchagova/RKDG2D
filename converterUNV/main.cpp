#include <iostream>
#include "FileConverter.h"

using namespace std;



// ----------------------------------------------------------------------------


int main()
{
    FileConverter converter("Mesh_Triangle.unv","Mesh2D");

    converter.importUNV();
    converter.exportRKDG();

    cout << "END" << endl;
    return 0;
}
