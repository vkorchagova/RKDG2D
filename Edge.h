#ifndef EDGE_H
#define EDGE_H

#include "numvector.h"
#include "RKDGCell.h"
#include <vector>

namespace std
{

class RKDGCell;

class Edge
{

public:

    //- Nodes edge consists of
    numvector<double, 2> nodes;

    //- Neighbour cell for edge
    numvector<RKDGCell*, 2> neibCells;

    //- Number of gauss points for edge
    int nGP = 2; //в глобальной переменной прописать веса и положения гауссовых точек

    //- Local numerical fluxes for edge
    numvector<double,2> localFluxes; // значения потоков в гауссовых точках - привязаны к ребру

public:
    Edge();
    ~Edge();

    //- Calculate local fluxes for edge
    getLocalFluxes();

}; // for Edge

} //for namespace std

#endif // EDGE_H
