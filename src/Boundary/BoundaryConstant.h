#ifndef BOUNDARYCONSTANT_H
#define BOUNDARYCONSTANT_H

#include "Boundary.h"
    
/// Fixed value boundary condition (constant)
///
/// Provides the constant value of conservative variables on the boundary.
/// 
/// Conservative variables are:
/// -   density;
/// -   normal component of velocity;
/// -   tangential component of velocity;
/// -   z-component of velocity (0);
/// -   volumetric total energy.

class BoundaryConstant : public Boundary
{
    /// Constant values of conservative variables on the boundary
    numvector<double, 5> fixedValue;

public:

    /// Constructor
    ///
    /// @param defValue values of conservative variables on the boundary
    BoundaryConstant(const numvector<double, 5>& defValue) : Boundary(), fixedValue(defValue) { type = "constant"; }

    /// Apply boundary conditions
    ///
    /// This function returns constant values of conservative variables on the boundary.
    ///
    /// @param  sol solution near the point n inside boundary cell (not used)
    /// @param  n   Normal direction related to boundary edge (not used)
    ///
    /// @return the fixed value of solution outside boundary cell 
    ///
    /// > By default the positive value of normal component of velocity directs outside flow region!
    /// > Here it should be directed into the flow domain (like an inlet)!
    /// > Therefore, normal component of velocity in this function is multiplied to -1.
    numvector<double, 5> applyBoundary(const numvector<double, 5>& solLeft = {0.0, 0.0, 0.0, 0.0, 0.0}, const Point& n = Point({0.0,0.0})) const override;

};

#endif // BOUNDARYCONSTANT_H
