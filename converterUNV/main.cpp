#include <iostream>
#include "FileConverter.h"

using namespace std;



// ----------------------------------------------------------------------------


int main()
{
    FileConverter converter("Mesh_1.unv","Mesh2D");

    converter.importUNV();

    cout << "END" << endl;
    return 0;
}
