#ifndef BUFFERS_H
#define BUFFERS_H

#include <map>
#include <vector>
#include "numvector.h"
#include "Params.h"

class Buffers
{
public:

    /// Buffers

    /// MPI buffer for full solution
    std::vector<double> forFullSOL; // size = nRealCellsTotal * dimS

    /// MPI buffers for local solution packs
    std::vector<double> forSendLocalSOL; // size = nRealCells * dimS
    std::vector<double> forRecvLocalSOL; // size = nRealCells * dimS

    /// MPI buffers for proc bound packs
    // std::map<int, std::vector<double>> forSendBoundSOL; // size = nProcPatches | nCellsOnBound * dimS
    std::vector<double> forSendBoundSOL; // size = prod(len(procPatches)) * dimS
    std::vector<double> forRecvBoundSOL; // size = nProcPatches | nCellsOnBound * dimS

    std::vector<double> forSolExport;
	std::vector<double> forSolExportRecv;


    /// Maps

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

