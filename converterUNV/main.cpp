#include <iostream>
#include "FileConverter.h"
#include <cstdlib>

using namespace std;


// ----------------------------------------------------------------------------


int main()
{

    FileConverter converter("square1_96.unv","mesh2D");


    converter.importUNV();
    converter.exportRKDG();
    
    cout << "END" << endl;
    return 0;
}
