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

#include "defs.h"		//- Basic arithmetics
//#include "Params.h"		//- All the manually defining parameters of the method
#include "Mesh.h"		//- The mesh
#include "Basis.h"		//- All about the basis functions for DG
#include "compService.h"//- Integration
#include "Physics.h"	//- The physical models
#include "Problem.h"	//- Initial-boundary staff
#include "Solution.h"   //- Solution storage
#include "Boundary.h"
#include "LimiterFinDiff.h"
#include "LimiterBJ.h"
#include "FluxLLF.h"		//- All about the flux evaluating
//#include "Indicator.h"
//#include "Limiter.h"	//- All about the monotonization
#include "Solver.h"		//- The whole spatial discretization module
#include "TimeControl.h"
#include "Writer.h"
#include "RungeKutta.h"

//#include "intel64\include\mpi.h"
#include <mpi.h>
#include <omp.h>

//- Proc rank 
int myRank;

//- Total number of procs
int numProcsTotal;

//- Status
MPI_Status status;

//- Request
MPI_Request request;


using namespace std;

int main(int argc, char* argv[])
{
    
    //int rank, size, ibeg, iend;
	//MPI_Status stat;

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &numProcsTotal);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    //cout << "size = " << numProcsTotal << "; rank = " << myRank << endl;

    //-----------------------

    string meshFileName = "mesh2D." + to_string(myRank);
    Buffers buf;
    Mesh mesh(meshFileName, buf);
    Basis basis(mesh.cells);
    Solution solution(basis);
    
    Physics physics;
    FluxLLF flux(physics);
    
    Writer writer(mesh, solution, physics);

    int order = 2;
    double tStart = 0.0;
    double tEnd = 1.0;
    double initTau = 2.5e-5;
    double outputInterval = 0.05;
    TimeControl time(mesh, tStart, tEnd, initTau, outputInterval);

    CaseInit caseName = SodCircle;
    Problem problem(caseName, mesh, time);

    Solver solver(basis, mesh, solution, problem, physics, flux, buf);

    LimiterBJ limiter(mesh.cells, solution);

    RungeKutta RK(order, basis, solver, solution, problem.bc, limiter, time);

    //-----------------------

    // create folder for solution
    #if !defined(__linux__)
        _mkdir("alphaCoeffs");
    #else
        mkdir("alphaCoeffs", S_IRWXU | S_IRGRP | S_IROTH);
    #endif

    // Set initial conditions
    if (tStart > 1e-10)
    {
        solver.setDefinedCoefficients("alphaCoeffs/" + to_string(tStart) + ".dat");
    }
    else
    {
    // get initial conditions
        solver.setInitialConditions();      // FIX GRAIENTS FOR nSHAPES > 3

        for (size_t i = 0; i < mesh.patches.size(); ++i)
            problem.bc[i]->applyBoundary(solution.SOL);

        solver.dataExchange();

        limiter.limit(solution.SOL);

        solver.collectSolution();
            
        if (myRank == 0) 
            writer.exportNativeCoeffs("alphaCoeffs/" + to_string(tStart) + ".dat");
    }

    physics.cpcv = problem.cpcv;
    //cout << physics.cpcv << endl;

    //for (int i = 0; i < mesh.patches.size(); ++i)
        //cout << "BC type #" << i << ":" << problem.bc[i]->type << endl;

    //for (size_t i = 0; i < mesh.cells.size(); ++i)
     //   cout << "sol #" << i << ": " << solution.SOL[i] << endl;

    //-----------------------

    
    //solver.assembleRHS(solution.SOL);
    
    double t0, t1;

    //for (double t = tStart; t < tEnd; t += tau)
    while (time.running())
    {
        if (myRank == 0) 
        {
            cout << "-------" << endl;
            cout << "Time = " << time.getTime() << " s" << endl;
            cout << "\ttau = " << time.getTau() << " s" << endl;
        }

        t0 = MPI_Wtime();
        RK.Tstep();
        t1 = MPI_Wtime();

        if (myRank == 0) 
        {
            cout << "\tstep time = " << t1 - t0 << " s" << endl;
        }

        if (time.isOutput())
        {
            solver.collectSolution();
            
            if (myRank == 0) 
                writer.exportNativeCoeffs("alphaCoeffs/" + to_string(time.getTime()) + ".dat"); // do it only on proc rank = 0
        }
    }

    /// separate program for visualisation: 
    /// load full mesh (prepare it in decomposer) and full pack of solution coeffs 
    /// and convert it to VTK
    

    MPI_Finalize();

    cout << "THE END" << endl;
	
	return 0;
}

//ibeg = rank*n + 1;
//iend = (rank + 1)*n;

//for (int i = ibeg; i <= ((iend > MAX) ? MAX : iend); ++i)
//{
	//	printf("process %d, %d^2=%d\n", rank.i.i * 1);
//}

/*if (rank == 1)
MPI_Send(&ibuf, 1, MPI_INT, 0, 5, MPI_COMM_WORLD);
if (rank == 0)
MPI_Recv(&ibuf, 1, MPI_INT, 1, 5, MPI_COMM_WORLD, &stat);*/

//if (rank == 0)
//	printf("process %d, TotSum =  %f\n", rank, res);
