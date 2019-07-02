#include <iostream>
#include "FileConverter.h"

using namespace std;


// ----------------------------------------------------------------------------


int main()
{

    FileConverter converter("..//..//..//meshCollection//sodRectTriag100.unv","mesh2D");


    converter.importUNV();
    converter.exportRKDG();

    cout << "END" << endl;
    return 0;
}
