#ifndef RKDGCELL_H
#define RKDGCELL_H

#include "numvector.h"
#include "Edge.h"
#include <vector>

namespace std
{

class Edge;

class RKDGCell
{
public:

    //- Edges cell consists of
    numvector<Edge*, 4> edges;

    //- Number of 2D Gauss points for cell
    int nGP = 4;

    //- Area of cell
    double Area;

    //- Mass center of cell
    double center;

    //- Calculate coordinates of cell nodes
    getCellCoordinates();

public:
    RKDGCell();
    ~RKDGCell();

}; // RKDGCell

} // End namespace std

#endif // RKDGCELL_H
