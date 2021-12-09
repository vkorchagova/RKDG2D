#ifndef PROBES_H
#define PROBES_H

#include "Mesh.h"
#include "Physics.h"
#include "Solution.h"
#include <fstream>

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

///
/// Probes of solution to files
///

class Probes
{
    /// Array of points for probes
    std::vector<Point> probes;

    /// Numbers of cells contains probe points
    std::vector<int> iCellsForProbes;

    /// Pressure values in probes
    std::vector<double> pForProbes;
    std::vector<double> pForProbesRecv;

private:

    /// Const reference to mesh
    const Mesh& mesh;

    /// Const reference to solution
    const Solution& solution;

    /// Const reference to physics
    const Physics& physics;

    /// File stream for 
    std::ofstream probesStream;

    /// Collect full solution 
    void collectFullSolution();

public:

    /// Constructor
    Probes(const Mesh& msh, const Solution& sln, const Physics& phs, const std::vector<Point>& _probes);

    /// Destructor
    ~Probes();

    //----- export
    void writeProbes(double t);

};


#endif // Probes_H