#ifndef EDGEBOUNDARYDIAGPROJECTION_H
#define EDGEBOUNDARYDIAGPROJECTION_H

#include "Point.h"
#include "EdgeBoundary.h"
#include "Mesh2D.h"

class Mesh2D;

class EdgeBoundaryDiagProjection : public EdgeBoundary
{
    //- Get projection of the defined point onto diagonal (only for uniform meshes!!!)
    Point getDiagProjection(const Point &p) const;

    //- Pointer to mesh
    const Mesh2D* mesh;

    //- Get number of cell in diagonal where given point placed
    int getNumDiagCell(const Point& p);

public:

    //- Default constructor
    EdgeBoundaryDiagProjection() : EdgeBoundary()  {}

    //- Construct using two nodes
    EdgeBoundaryDiagProjection(const Node& p1, const Node& p2) : EdgeBoundary(p1, p2) {}

    //- Destructor
    virtual ~EdgeBoundaryDiagProjection() = default;

    //// RKDG methods

    //- Set BC function
    virtual void setBoundary(const std::shared_ptr<Boundary>& bound) override {}

    //- Set mesh pointer
    void setMeshPointer(const Mesh2D& msh) { mesh = &msh; }

    //- Calculate local fluxes in gauss points
    virtual void getLocalFluxes(const Flux& flux) override;

    //- Compute max speed on edge
    virtual void getMaxUL() override {}

};

#endif // EDGEBOUNDARYDIAGPROJECTION_H
