#ifndef INDICATOR_H
#define INDICATOR_H

#include <vector>
#include "Mesh2D.h"
#include "numvector.h"


class Indicator
{
public:

    //- Mesh
    const Mesh2D& mesh;

    //- Problem
    const Problem& problem;
    
    //- List of troubed cells
    std::vector<int> tCells;
    
public:

    //- Constructor
    Indicator(const Mesh2D& msh, const Problem& prb);

    //- Destructor
    virtual ~Indicator() {}
    
    //- Check discontinuities
    virtual std::vector<int> checkDiscontinuities() const = 0;
    
    //- Write troubled cells in VTK format
    //void writeTroubledCellsVTK();
    
};




#endif
