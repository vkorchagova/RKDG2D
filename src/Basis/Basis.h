#ifndef BASIS_H
#define BASIS_H

#include <vector>
#include <function>
#include "Point.h"
///#include <Mesh.h>

/// Form functions for mesh cells

class Basis
{
    //- Constant reference to mesh
    ///const Mesh& mesh;
public:

    Basis() {};
    //Basis(int nF) { nShapes = nF;};
    
    //- Coefficients for form functions
    vector<numvector<double, nShapes>> phiCoeffs;
    
    //- List of form functions as is
    std::vector<std::function<double(const Point&)>> phi;
    
    //- Gramian matrix
    numvector< numvector<double,nShapes>, nShapes> gramian;

};

#endif // BASIS_H

