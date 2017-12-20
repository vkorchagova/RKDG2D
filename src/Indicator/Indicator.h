#ifndef INDICATOR_H
#define INDICATOR_H

#include <vector>
#include "Mesh2D.h"
#include "numvector.h"


class Indicator
{
protected:

    //- Mesh
    const Mesh2D& mesh;
    
    //- Coeffs
    const std::vector<numvector<double, 5 * nShapes>> alpha;
    
public:

    //- Constructor
    Indicator(const Mesh2D& msh, const std::vector<numvector<double, 5 * nShapes>>& coeffs) : mesh(msh), alpha(coeffs) {};
    
    //- Check discontinuities
    virtual std::vector<double> checkDiscontinuities() const = 0;

};




#endif