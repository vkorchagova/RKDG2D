#ifndef INDICATOR_H
#define INDICATOR_H

#include <vector>
#include "Mesh2D.h"
#include "numvector.h"

///
/// Abstract class for troubled cell indicators
///

class Indicator
{
public:

    /// Constant reference to mesh
    const Mesh2D& mesh;

    /// Constant reference to problem
    const Problem& problem;
    
    /// List of troubled cells
    std::vector<int> tCells;
    
public:

    /// Constructor
    Indicator(const Mesh2D& msh, const Problem& prb);

    /// Destructor
    virtual ~Indicator() {}
    
    /// Check discontinuities
    virtual std::vector<int> checkDiscontinuities() const = 0;
    
};

#endif // INDICATOR_H
