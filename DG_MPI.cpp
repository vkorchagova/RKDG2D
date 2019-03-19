// DG_MPI.cpp
//


#if !defined(__linux__)
#include <direct.h>
#endif

#include <sys/stat.h>

#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>
#include <time.h>
#include "numvector.h"
#include "Buffers.h"

#include "defs.h"		/// Basic arithmetics
//#include "Params.h"		/// All the manually defining parameters of the method
#include "Mesh.h"		/// The mesh
#include "Basis.h"		/// All about the basis functions for DG
#include "compService.h"/// Integration
#include "Physics.h"	/// The physical models
#include "Problem.h"	/// Initial-boundary staff
#include "Solution.h"   /// Solution storage
#include "Boundary.h"
#include "LimiterFinDiff.h"
#include "LimiterBJ.h"
#include "LimiterWENOS.h"
#include "FluxLLF.h"		/// All about the flux evaluating
#include "FluxHLL.h"
#include "FluxHLLC.h"
//#include "Indicator.h"
//#include "Limiter.h"	/// All about the monotonization
#include "Solver.h"		/// The whole spatial discretization module
#include "TimeControl.h"
#include "Writer.h"
#include "RungeKutta.h"

//#include "intel64\include\mpi.h"
#include <mpi.h>
#include <omp.h>

/// Proc rank 
int myRank;

/// Total number of procs
int numProcsTotal;

/// Status
MPI_Status status;

/// Debug
bool debug;

/// Log file to save data
std::ofstream logger;


using namespace std;

int main(int argc, char* argv[])
{
    
    //int rank, size, ibeg, iend;
    //MPI_Status stat;


    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &numProcsTotal);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    debug = (true && myRank == 0); // if you want save log type true else type false
    logger.open("timeStat." + to_string(numProcsTotal));

    //if (myRank == 0) cout << "size = " << numProcsTotal << "; rank = " << myRank << endl;

    // create folder for solution
    #if !defined(__linux__)
        _mkdir("alphaCoeffs");
    #else
        mkdir("alphaCoeffs", S_IRWXU | S_IRGRP | S_IROTH);
    #endif

    ///----------------------

    CaseInit caseName = SodCircle;

    double tStart = 0.0;
    double tEnd = 0.001;
    double initTau = 1e-3;
    double outputInterval = 0.2;

    int order = 2;


    string meshFileName = "mesh2D"; 
    if (numProcsTotal > 1)
        meshFileName += "." + to_string(myRank);

    Buffers buf;

    Mesh mesh(meshFileName, buf);
    Mesh fullMesh(buf);

    if (myRank == 0)
    {
        meshFileName = "mesh2D";
        fullMesh.importMesh(meshFileName);
    }

    Basis basis(mesh.cells);
    Solution solution(basis);
    
    Physics physics;
    FluxHLL flux(physics);
    
    Writer writer(fullMesh, solution, physics);
    TimeControl time(mesh, tStart, tEnd, initTau, outputInterval);

    Problem problem(caseName, mesh, time);
    Solver solver(basis, mesh, solution, problem, physics, flux, buf);

    LimiterWENOS limiter(mesh.cells, solution, physics);
    RungeKutta RK(order, basis, solver, solution, problem.bc, limiter, time);

    ///----------------------

    //writer.exportMeshVTK("mesh2D.vtk");

    // Set initial conditions
    if (tStart > 1e-10)
    {
        solver.setDefinedCoefficients("alphaCoeffs/" + to_string(tStart) + ".dat");
    }
    else
    {
    // get initial conditions
        solver.setInitialConditions();      // FIX GRAIENTS FOR nSHAPES > 3

        // for (size_t i = 0; i < mesh.patches.size(); ++i)
        //     problem.bc[i]->applyBoundary(solution.SOL);

        solver.dataExchange();

        limiter.limit(solution.SOL);

        solver.collectSolution();
        solver.collectSolutionForExport();
            
        if (myRank == 0) 
        {
            writer.exportNativeCoeffs("alphaCoeffs/" + to_string(tStart) + ".dat");
            writer.exportFrameVTK("alphaCoeffs/" + to_string(tStart) + ".vtk");
        }
    }

    physics.cpcv = problem.cpcv;
    //cout << physics.cpcv << endl;

    //for (int i = 0; i < mesh.patches.size(); ++i)
        //cout << "BC type #" << i << ":" << problem.bc[i]->type << endl;

    //for (size_t i = 0; i < mesh.cells.size(); ++i)
     //   cout << "sol #" << i << ": " << solution.SOL[i] << endl;

    ///----------------------
    
    double t0, t1;

    double meanCpuTime = 0.0;
    double totalCpuTime = 0.0;
    int nSteps = 0;
    

	/// THE MAIN CYCLE THROUGH THE TIME!
    while (time.running())
    {
        if (myRank == 0) 
        {
            cout << "-------" << endl;
            cout << "Time = " << time.getTime() << " s" << "\t\t";
            cout << "Tau = " << time.getTau() << " s" << "\t\t";
        }

        t0 = MPI_Wtime();
        RK.Tstep();
        t1 = MPI_Wtime();

        if (debug) logger << "RK.Tsep(): " << t1 - t0  << "\n----------" << endl;

        if (myRank == 0) 
        {
            cout << "CPU time = " << t1 - t0 << " s" << endl;
            totalCpuTime += t1 - t0;
            nSteps ++;
        }

        if (time.isOutput())
        {
            solver.collectSolution();
            solver.collectSolutionForExport();
            
            if (myRank == 0) 
            {
                writer.exportNativeCoeffs("alphaCoeffs/" + to_string(time.getTime()) + ".dat");
                writer.exportFrameVTK("alphaCoeffs/" + to_string(time.getTime()) + ".vtk");
            }
        }
    }

    if (myRank == 0)
    {
        cout << "============" << endl;
        cout << "Total CPU time: " << totalCpuTime << " s" << endl;
        cout << "Number of time steps: " << nSteps << endl;
        cout << "Mean CPU timestep: " << totalCpuTime / double(nSteps) << " s" << endl;
        cout << "============" << endl;
        cout << "THE END" << endl;
    }

    logger.close();


    MPI_Finalize();
	
	return 0;
}

