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
#include "Mesh.h"		//- The very mesh
//#include "Basis.h"		//- All about the basis functions for DG
#include "compService.h"//- Integration
//#include "Physics.h"	//- The physical models
//#include "Problem.h"	//- Initial-boundary staff
//#include "Solution.h"   //- Solution storage
//#include "Boundary.h"
//#include "Flux.h"		//- All about the flux evaluating
//#include "Indicator.h"
//#include "Limiter.h"	//- All about the monotonization
//#include "Solver.h"		//- The whole spatial discretization module
//#include "TimeControl.h"
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
