#ifndef BOUNDARYNOZZLEINLET_H
#define BOUNDARYNOZZLEINLET_H

#include "Boundary.h"
#include "Problem.h"

class Problem;

class BoundaryNozzleInlet : public Boundary
{
    //- Reference to problem
    const Problem& problem;
    
    //- Total pressure, Pa
    double pTot;
    
    //- Temperature (fixed), K
    double T;
    
    //- Molar mass, kg / mol (input: g/mol)
    double M;
    
    //- Universal gas constant, J / mol / K)
    double R0;

public:
    //BoundaryNozzleInlet() : Boundary() { type = "slip";}
    
     //- Constructor
    BoundaryNozzleInlet(double p0, double T, double M, const Problem& prb, const numvector<double,5>& u0 = {1,0,0,0,1}) 
    : Boundary(), pTot(p0), T(T), M(M), problem(prb) { type = "nozzle_inlet"; R0 = 10e3*8.314459848;}
    ~BoundaryNozzleInlet() {};
    
    //- Apply boundary condition
    numvector<double, 5> applyBoundary(const numvector<double, 5>& solLeft = {0.0, 0.0, 0.0, 0.0, 0.0}, const Point& n = Point({0.0,0.0}), int numGP = 0) const override;
};

#endif // BOUNDARYNOZZLEINLET_H
