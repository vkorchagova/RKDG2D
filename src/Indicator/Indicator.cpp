#include "Indicator.h"
#include <fstream>

using namespace std;
                                                
Indicator::Indicator (const Mesh& msh, const Solution& sln) : mesh(msh), solution(sln)
{    
    values.resize(mesh.nRealCells*dimPh);
}