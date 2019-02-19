#ifndef WRITER_H
#define WRITER_H

#include "Mesh.h"
#include "Physics.h"
#include "Solution.h"

// TODO: move here mesh. ExportVtkPoly & solution.export as exportNativeCoeffs

/// Proc rank 
extern int myRank;

/// Size of ...
extern int numProcsTotal;

/// Status
extern MPI_Status status;

/// Debug
extern bool debug;

/// Log file to save data
extern std::ofstream logger;


class Writer
{
    /// Const reference to mesh
    const Mesh& mesh;

    /// Const reference to solution
    const Solution& solution;

    /// Const reference to physics
    const Physics& physics;

public:

    /// Constructor
    Writer(const Mesh& msh, const Solution& sln, const Physics& phs) : mesh(msh), solution(sln), physics(phs) {}

    /// Export mesh to VTK
    void exportMeshVTK(const std::string& fileName) const;
    void exportMeshVTK(std::ostream& wStream) const;

    /// Export solution to VTK
    void exportFrameVTK(const std::string& fileName) const;
    void exportFrameVTK(std::ostream& wStream) const;

    /// Export solution coeffs 
    void exportNativeCoeffs(const std::string& fileName) const;
    void outputNativeCoeffs(std::ostream& wStream) const;
    void outputFullNativeCoeffs(std::ostream& wStream) const;
    void outputFullNativeCoeffs(std::ostream& wStream, std::vector<numvector<double, dimS>> coeffs) const;

    /// MPI

    /// Collect full solution 
    void collectFullSolution();

};


#endif // WRITER_H