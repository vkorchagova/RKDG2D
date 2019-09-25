#ifndef INDICATOR_H
#define INDICATOR_H

#include <vector>
#include "Mesh.h"
#include "Solution.h"
#include "numvector.h"

///
/// Abstract class for indicator of problem cells
///

class Indicator
{
 
protected:

    /// Constant reference to mesh
    const Mesh& mesh;

    /// Constant reference to solution
    const Solution& solution;
    
public:

    /// Constructor
    Indicator (const Mesh& msh, const Solution& sln) : mesh(msh), solution(sln) {}
    
    /// Destructor
    virtual ~Indicator() {}
    
    /// Check discontinuities
    virtual std::vector<int> checkDiscontinuities() const = 0;
    
    /// Write troubled cells in VTK format
    //void writeTroubledCellsVTK();
    
};




#endif
