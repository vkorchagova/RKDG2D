// DG_MPI.cpp: main
//

//#include "stdafx.h"
//#include <stdio.h>
//#include <iostream>
//#include <iostream>
//#include <iomanip>
//#include <fstream>
//#include <cmath>
//#include <algorithm>
//#include <vector>
//#include <string>
//#include <time.h>
//#include "numvector.h"

//#include "defs.h"
//#include "Params.h"

//#include "Mesh.h"
//#include "NEX.h"
//#include "Physics.h"
//#include "Problem.h"
//#include "Boundary.h"
//#include "Flux.h"
//#include "Indicator.h"
//#include "Limiter.h"
//#include "Solver.h"
//#include "TimeControl.h"
//#include "TimeStepper.h"

#include <mpi.h>
//#include "intel64\include\mpi.h"
//#include <omp.h>

using namespace std;

int main(int argc, char* argv[])
{
    
	int rank, size, ibeg, iend;
	MPI_Status stat;

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	MPI_Finalize();
	
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