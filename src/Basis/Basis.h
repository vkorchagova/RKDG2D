#ifndef BASIS_H
#define BASIS_H

#include <vector>
#include <functional>
#include "Point.h"

#include "Mesh.h"
#include "Params.h"



/// Form functions for mesh cells

class Basis
{
    /// Constant reference to mesh
    const std::vector<std::shared_ptr<Cell>>& cells;

    /// Init basis functions
    void initBasisFunctions();

    /// Init gramian matrices
    void initGramian();

public:
    //Basis(int nF) { nShapes = nF;};

    Basis(const std::vector<std::shared_ptr<Cell>>& cells);
    
    /// Coefficients for form functions
    std::vector<numvector<double, nShapes>> phiCoeffs;
    
    /// List of form functions as is
    std::vector<std::function<double(int iCell, const Point&)>> phi;
	//std::vector<std::function<double(const Point&)>> phi;

    /// List of form functions' gradients
    std::vector<std::function<Point(int iCell, const Point&)>> gradPhi;
    
    /// Gramian matrix
    /// Typical look for my basis:
    /// 1     0         0
    /// 0 phi1*phi1 phi1*phi2
    /// 0 phi2*phi1 phi2*phi2
    ///
    /// The good idea is to save only one triangle of symmetric matrix.
    /// NOTE: in case of constant functions will be empty.
    /// NOTE: gramians for each cells are not similar in unstructured mesh!
    ///
    /// Gramian form for one cell in case of linear functions: {g[1][1], g[1][2], g[2][2]}
	std::vector< std::vector<std::vector<double>> > gramian;

    /// Get coefficients of projection of function foo onto cell basis
    numvector<double, dimS> projection(const std::function<numvector<double, dimPh>(const Point& point)>& init, int iCell) const;

    /// Correct alpha coeffs in case of non-orthogonal basis functions
    numvector<double, dimS> correctNonOrthoCell(const numvector<double, dimS>& rhs, const std::vector<std::vector<double>>& g) const;

};

#endif // BASIS_H

