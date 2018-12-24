// DG_MPI.cpp
//


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
//#include "Boundary.h"
#include "FluxLLF.h"		//- All about the flux evaluating
//#include "Indicator.h"
//#include "Limiter.h"	//- All about the monotonization
#include "Solver.h"		//- The whole spatial discretization module
#include "TimeControl.h"
#include "Writer.h"
//#include "TimeStepper.h"

//#include "intel64\include\mpi.h"
#include <mpi.h>


using namespace std;

int main(int argc, char* argv[])
{
    
    int rank, size, ibeg, iend;
	MPI_Status stat;

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    cout << "size = " << size << "; rank = " << rank << endl;

    string meshFileName = "mesh2D";// + to_string(rank);

    //-----------------------

    Mesh mesh(meshFileName);
    Basis basis(mesh.cells);
    TimeControl time(mesh);
    Physics physics;
    Problem problem(SodX,mesh,time);
    Solution solution(basis);
    FluxLLF flux(physics);
    Solver solver(basis, mesh, solution, problem, flux);
    Writer writer(mesh, solution, physics);

    solver.setInitialConditions();

    physics.cpcv = problem.cpcv;
    cout << physics.cpcv << endl;

    for (int i = 0; i < mesh.patches.size(); ++i)
        cout << "BC type #" << i << ":" << problem.bc[i]->type << endl;

    //-----------------------

    writer.exportMeshVTK("mesh2D.vtk");
    writer.exportFrameVTK("0.vtk");
    writer.exportNativeCoeffs("0.dat");

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
