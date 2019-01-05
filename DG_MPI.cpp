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


using namespace std;

int main(int argc, char* argv[])
{
    
    int rank, size, ibeg, iend;
	MPI_Status stat;

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    cout << "size = " << size << "; rank = " << rank << endl;

    //-----------------------

    string meshFileName = "mesh2D";// + to_string(rank);
    Mesh mesh(meshFileName);
    Basis basis(mesh.cells);
    Solution solution(basis);
    
    Physics physics;
    FluxLLF flux(physics);
    
    Writer writer(mesh, solution, physics);

    int order = 2;
    double tStart = 0.0;
    double tEnd = 0.2;
    double initTau = 1e-4;
    TimeControl time(mesh, tStart, tEnd, initTau);

    CaseInit caseName = SodCircle;
    Problem problem(caseName, mesh, time);

    Solver solver(basis, mesh, solution, problem, physics, flux);

    LimiterFinDiff limiter;

    RungeKutta RK(order, basis, solver, solution, problem.bc, limiter, time);

    //-----------------------

    // create folder for solution
    #if !defined(__linux__)
        _mkdir("alphaCoeffs");
    #else
        mkdir("alphaCoeffs", S_IRWXU | S_IRGRP | S_IROTH);
    #endif

    // get initial conditions
    solver.setInitialConditions();      // FIX GRAIENTS FOR nSHAPES > 3
    writer.exportFrameVTK("alphaCoeffs/0.vtk");
    writer.exportNativeCoeffs("alphaCoeffs/0.dat");

    physics.cpcv = problem.cpcv;
    cout << physics.cpcv << endl;

    for (int i = 0; i < mesh.patches.size(); ++i)
        cout << "BC type #" << i << ":" << problem.bc[i]->type << endl;

    //for (size_t i = 0; i < mesh.cells.size(); ++i)
     //   cout << "sol #" << i << ": " << solution.SOL[i] << endl;

    //-----------------------

    // 1st step: apply boundary
    for (size_t i = 0; i < mesh.patches.size(); ++i)
    {
        problem.bc[i]->applyBoundary(solution.SOL);
    }

    limiter.limit(solution.SOL);
    //solver.assembleRHS(solution.SOL);
    
    solution.SOL_aux=solution.SOL;

    
    double t0, t1;

    //for (double t = tStart; t < tEnd; t += tau)
    while (time.running())
    {
        cout << "-------" << endl;
        cout << "Time: " << time.getTime() << endl;
        cout << "Time step: " << time.getTau() << endl;

        t0 = MPI_Wtime();
        RK.Tstep();
        t1 = MPI_Wtime();

        cout << "\tstep time = " << t1 - t0 << " s" << endl;


        writer.exportFrameVTK("alphaCoeffs/" + to_string(time.getTime()) + ".vtk"); 
        writer.exportNativeCoeffs("alphaCoeffs/" + to_string(time.getTime()) + ".dat");
    }

    
    

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
