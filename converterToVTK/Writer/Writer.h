#ifndef WRITER_H
#define WRITER_H

#include "Mesh.h"
#include "Physics.h"
#include "Solution.h"

// TODO: move here mesh. ExportVtkPoly & solution.export as exportNativeCoeffs


class Writer
{
    //- Const reference to mesh
    const Mesh& mesh;

    //- Const reference to solution
    const Solution& solution;

    //- Const reference to physics
    const Physics& physics;

public:

    //- Constructor
    Writer(const Mesh& msh, const Solution& sln, const Physics& phs) : mesh(msh), solution(sln), physics(phs) {}

    //- Export mesh to VTK
    void exportMeshVTK(const std::string& fileName) const;
    void exportMeshVTK(std::ostream& wStream) const;

    //- Export solution to VTK
    void exportFrameVTK(const std::string& fileName) const;
    void exportFrameVTK(std::ostream& wStream) const;

};


#endif // WRITER_H