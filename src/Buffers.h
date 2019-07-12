#ifndef BUFFERS_H
#define BUFFERS_H

#include <map>
#include <vector>
#include "numvector.h"
#include "Params.h"

///
/// Buffers and maps for MPI operations
///

class Buffers
{
public:

    //----- Buffers

    /// MPI buffer for full solution
    /// size = nRealCellsTotal * dimS
    std::vector<double> forFullSOL; 

    /// MPI buffer for sending of local solution pack in internal cells
    /// size = nRealCells * dimS
    std::vector<double> forSendLocalSOL;

    /// MPI buffer for receiving of local solution pack in internal cells
    /// size = nRealCells * dimS
    std::vector<double> forRecvLocalSOL;

    /// MPI buffer for sending of proc bound pack
    /// size = prod(len(procPatches)) * dimS
    std::vector<double> forSendBoundSOL;

    /// MPI buffer for receiving of proc bound pack
    /// size = prod(len(procPatches)) * dimS 
    std::vector<double> forRecvBoundSOL;

    /// MPI buffer for sending of local pack of solution in cell centres
    std::vector<double> forSolExport;

    /// MPI buffer for receiving of local pack of solution in cell centres
	std::vector<double> forSolExportRecv;

    //----- Maps

    /// Number of "real" cells on each process
    std::vector<int> nCellsPerProc;

    /// Number of solution coeffs on each proc
    std::vector<int> nCoeffsPerProc;
    std::vector<int> nSolPerProcExp; 
    std::vector<int> boundProcSolSizesS;
    std::vector<int> boundProcSolSizesR;

    /// Array of displacements for MPI_Gatherv
    std::vector<int> mpiDisplMesh; 
    std::vector<int> mpiDisplSOL;    
    std::vector<int> mpiDisplSolExp; 
    std::vector<int> boundProcSolDisplS;
    std::vector<int> boundProcSolDisplR;

    /// Full global numeration of cells
    std::vector<int> globalMap;

    /// Default constructor
    Buffers()  {};

    /// Destructor
    ~Buffers() {};

};

#endif

