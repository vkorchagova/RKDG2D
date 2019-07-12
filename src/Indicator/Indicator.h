#ifndef INDICATOR_H
#define INDICATOR_H

#include <vector>
#include "Mesh.h"
#include "numvector.h"

///
/// Abstract class for indicator of problem cells
///

class Indicator
{
public:

    /// Constant reference to mesh
    const Mesh& mesh;

    /// Problem
    //const Problem& problem;
    
    /// List of troubed cells
    //std::vector<int> tCells;
    
public:

    /// Constructor
    Indicator (const Mesh& msh) : mesh(msh) {}
    
    /// Destructor
    virtual ~Indicator() {}
    
    /// Check discontinuities
    virtual std::vector<int> checkDiscontinuities() const = 0;
    
    /// Write troubled cells in VTK format
    //void writeTroubledCellsVTK();
    
};




#endif
