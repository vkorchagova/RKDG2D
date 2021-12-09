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
#include "LimiterBJVertex.h"
#include "LimiterWENOS.h"
#include "LimiterRiemannWENOS.h"
#include "LimiterRiemannBJ.h"
#include "FluxLLF.h"		/// All about the flux evaluating
#include "FluxHLL.h"
#include "FluxHLLC.h"
#include "IndicatorEverywhere.h"
#include "IndicatorNowhere.h"
#include "IndicatorBJ.h"
#include "IndicatorShu.h"
//#include "Limiter.h"	/// All about the monotonization
#include "Solver.h"		/// The whole spatial discretization module
#include "SolverHybrid.h"
#include "TimeControl.h"
#include "Writer.h"
#include "RungeKutta.h"

#include "probes.h"

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

    debug = false;//(true && myRank == 0); // if you want save log type true else type false
    logger.open("timeStat." + to_string(numProcsTotal));

    //if (myRank == 0) cout << "size = " << numProcsTotal << "; rank = " << myRank << endl;

    // create folder for solution
    #if !defined(__linux__)
        _mkdir("alphaCoeffs");
    #else
        mkdir("alphaCoeffs", S_IRWXU | S_IRGRP | S_IROTH);
    #endif
	
    ///----------------------

    CaseInit caseName = Munday;//ForwardStep;//DoubleMach;//Munday;//ConvDivNozzle;

    double tStart = 0.0;
    double tEnd = 0.1;

    double initTau = 2e-7;
    double outputInterval = 1e-3;//0.5;//5e-2;

    bool isDynamic = true;
    double maxCo = 0.25;
    double maxTau = 1.0;
    double maxTauGrowth = 0.1; 

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
    // FluxHLLC flux(physics);

    FluxHLL flux1(physics);
    FluxHLLC flux2(physics);
    
    Writer writer(fullMesh, solution, physics);
    TimeControl time(mesh, physics, solution, tStart, tEnd, initTau, outputInterval, isDynamic, maxCo, maxTau, maxTauGrowth);
    Problem problem(caseName, mesh, time, physics);
    // Solver solver(basis, mesh, solution, problem, physics, flux, buf);
    SolverHybrid solver(basis, mesh, solution, problem, physics, flux1, flux2, buf);

    // IndicatorEverywhere indicator(mesh, solution);
    IndicatorBJ indicator(mesh, solution);
    LimiterRiemannWENOS limiter(mesh, solution, physics, indicator, buf);
    // LimiterWENOS limiter(mesh, solution, physics, indicator, buf);
    // LimiterBJ limiter(mesh, solution, physics, indicator, buf);
    // LimiterRiemannBJ limiter(mesh, solution, physics, indicator, buf);
    RungeKutta RK(order, basis, solver, solution, limiter, time);

    ///---------------------

    //writer.exportMeshVTK("mesh2D.vtk");

    // Set initial conditions
    if (tStart > 1e-10)
    {
        solver.restart("alphaCoeffs/" + to_string(tStart) + ".dat");
    }
    else
    {
    // get initial conditions
        solver.setInitialConditions();      // FIX GRAIENTS FOR nSHAPES > 3

        // for (size_t i = 0; i < mesh.patches.size(); ++i)
        //     problem.bc[i]->applyBoundary(solution.SOL);

        solver.dataExchange();

        limiter.limitSolution();

        if (myRank == 0) std::cout << "Initial conditions OK" << std::endl;

        solver.collectSolution();
        solver.collectSolutionForExport();
        solver.collectMeanSolutionForExport();
            
        if (myRank == 0) 
        {
            writer.exportNativeCoeffs("alphaCoeffs/" + to_string(tStart) + ".dat");
            writer.exportFrameVTK("alphaCoeffs/" + to_string(tStart) + ".vtk");
        }
    }
    
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

    vector<Point> microphones = 
    {
        Point({-0.08,0.1}),
        Point({-0.01,0.1}),
        Point({ 0.01,0.1}),
        Point({ 0.08,0.1}),
        Point({-0.08,0.2}),
        Point({-0.01,0.2}),
        Point({ 0.01,0.2}),
        Point({ 0.08,0.2})
    };

    Probes probeManager(mesh, solution, physics, microphones);
    

	/// THE MAIN CYCLE THROUGH THE TIME!
    while (time.running())
    {
        // make new time step

        t0 = MPI_Wtime();
        RK.Tstep();
        solver.computeTimeAvgSolution(time.getTau());
        t1 = MPI_Wtime();

        if (myRank == 0) 
        {
            cout << "-------" << endl;
            cout << "Time = " << time.getTime() << " s" << "\t\t";
            cout << "Tau = " << time.getTau() << " s" << "\t\t";
        }

        if (debug)
        {
            logger << "-------" << endl;
            logger << "Time = " << time.getTime() << " s" << "\t\t";
            logger << "Tau = " << time.getTau() << " s" << "\t\t";
        }

        if (debug) logger << "RK.Tstep(): " << t1 - t0  << "\n----------" << endl;

        if (myRank == 0) 
        {
            cout << "CPU dt time = " << t1 - t0 << " s" << "\t\t";
            totalCpuTime += t1 - t0;
            cout << "Execution time = " << totalCpuTime << " s" << endl;
            nSteps ++;
        }

        // check energy conservation

        double totalEnergy = 0.0;

//#pragma omp parallel for reduction(+:totalEnergy)
        //for (const shared_ptr<Cell> cell : mesh.cells)
        for ( int i = 0; i < mesh.nRealCells; ++i )
        {
            const shared_ptr<Cell> cell = mesh.cells[i];
            function<double(const Point&)> eTotal = [&](const Point& x)
            {
                return solution.reconstruct(cell->number, x, E);
            };

            totalEnergy += integrate(*(cell), eTotal);
        }
        
        if (numProcsTotal > 1)
        {
            double eTotal = 0.0;
            MPI_Reduce(
                &totalEnergy,
                &eTotal,
                1,
                MPI_DOUBLE,
                MPI_SUM,
                0,
                MPI_COMM_WORLD
            );
            
            //cout << myRank << ' ' << totalEnergy << endl;
            if (myRank == 0)
            {
                cout.precision(16);
                cout << "eTotal = " << eTotal << endl;
                cout.precision(6);
            }
        }
    	else
    	{
            if (myRank == 0)
            {
                cout.precision(16);
                cout << "eTotal = " << totalEnergy << endl;
                cout.precision(6);
            }
        }

        

        // write results in case of output time

        if (time.isOutput())
        {
			solver.collectSolution();
            solver.collectSolutionForExport();
            solver.collectMeanSolutionForExport();
            
            if (myRank == 0)
            {
                writer.exportNativeCoeffs("alphaCoeffs/" + to_string(time.getTime()) + ".dat");
                writer.exportFrameVTK("alphaCoeffs/" + to_string(time.getTime()) + ".vtk", time.getTime());
            }
        }

        if (totalEnergy != totalEnergy) //nan 
            exit(1);

        // update time step
        time.updateTimeStep();

        probeManager.writeProbes(time.getTime());
    }
    // cout << "Save last time point..." << endl;
    // // just for last time point
    // solver.collectSolution(); 
    // solver.collectSolutionForExport();
    
    if (myRank == 0)
    {
        writer.exportNativeCoeffs("alphaCoeffs/" + to_string(tEnd) + ".dat");
        writer.exportFrameVTK("alphaCoeffs/" + to_string(tEnd) + ".vtk");
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

    if (debug)
    {
        logger << "============" << endl;
        logger << "Total CPU time: " << totalCpuTime << " s" << endl;
        logger << "Number of time steps: " << nSteps << endl;
        logger << "Mean CPU timestep: " << totalCpuTime / double(nSteps) << " s" << endl;
        logger << "============" << endl;
        logger << "THE END" << endl;
    }

    logger.close();


    MPI_Finalize();
	
	return 0;
}

