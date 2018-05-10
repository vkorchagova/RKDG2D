//- RKDG 2D v.0.1
//  Structured rectangular mesh

#if !defined(__linux__)
#include <direct.h>
#endif

#include <sys/stat.h>

#include <stdio.h>
#include <iostream>
#include <string>

#include "defs.h"
#include "TimeClass.h"
#include "TimeControl.h"
#include <time.h>
#include "Mesh2D.h"
#include "Solver.h"
#include "FluxLLF.h"
#include "FluxHLL.h"
#include "FluxHLLC.h"
#include "IndicatorNowhere.h"
#include "IndicatorEverywhere.h"
#include "IndicatorKXRCF.h"
#include "IndicatorHarten.h"
#include "LimiterFinDiff.h"
#include "LimiterMUSCL.h"
#include "LimiterWENOS.h"


using namespace std;


int main(int argc, char** argv)
{    
    // Time parameters


    double Co = 0.1;
    double tEnd = 0.2;

    double initDeltaT = 1e-4;
    double maxDeltaT = 1.0;
    double maxTauGrowth = 1.2;
    bool isDynamicTimeStep = true;

    int freqWrite = 1;


    // ---------------

    // Initialize time
    Time time;

    // Initialize problem
    Problem problem(time);

    // Initialize mesh
    Mesh2D mesh("../RKDG2D/mesh2D",problem);

    // Set BC
    problem.setBoundaryConditions(mesh.patches);

    // Initialize flux
    FluxLLF numFlux(problem);

    // Initialize solver
    Solver solver(mesh, problem, numFlux);

    // Initialize time controller
    TimeControl dynamicTimeController(mesh,Co,maxDeltaT,maxTauGrowth,initDeltaT,isDynamicTimeStep);

    // Initialize indicator
    IndicatorKXRCF indicator(mesh, problem);

    // Initialize limiter
    LimiterWENOS limiter(indicator, problem);

    // ---------------

    #if !defined(__linux__)
        _mkdir("alphaCoeffs");
    #else
        mkdir("alphaCoeffs", S_IRWXU | S_IRGRP | S_IROTH);
    #endif

    cout << "---------\nt = " << 0 << endl;

    // Set initial conditions
    solver.setInitialConditions();

    limiter.limit(solver.alphaPrev);

    ofstream output;
    output.open("alphaCoeffs/0.000000.vtk");

    mesh.exportMeshVTK(output);

    output << "CELL_DATA " << mesh.nCells << endl;

    output << "FIELD rho 1" << endl;
    output << "rho" << " 1" << " " << mesh.nCells << " float" << endl;

    for (int i = 0; i < mesh.nCells; ++i)
       output << mesh.cells[i]->reconstructSolution(mesh.cells[i]->getCellCenter(),0) << endl;

    output << "FIELD rhoU 1" << endl;
    output << "rhoU " << " 3" << " " << mesh.nCells << " float" << endl;

    for (int i = 0; i < mesh.nCells; ++i)
    {
       output << mesh.cells[i]->reconstructSolution(mesh.cells[i]->getCellCenter(),1) << ' ';
       output << mesh.cells[i]->reconstructSolution(mesh.cells[i]->getCellCenter(),2) << ' ';
       output << 0.0 << endl;
    }

    output << "FIELD e 1" << endl;
    output << "e" << " 1" << " " << mesh.nCells << " float" << endl;

    for (int i = 0; i < mesh.nCells; ++i)
       output << mesh.cells[i]->reconstructSolution(mesh.cells[i]->getCellCenter(),4) << endl;

    output << "POINT_DATA " << mesh.nodes.size() << endl;

    output.close();

   // time step

    double tau = initDeltaT;

    vector<numvector<double, 5*nShapes>> lhs = solver.alphaPrev; // coeffs


    // run Runge --- Kutta 2 TVD

    vector<numvector<double, 5*nShapes>> k1, k2;

    k1.resize(mesh.nCells);
    k2.resize(mesh.nCells);

    clock_t t1, t2, t00;

    t00 = clock();

    int iT = 1; //iteration number

    for (double t = tau; t <= tEnd + 0.5*tau; t += tau)
    {
       t1 = clock();

       time.updateTime(t);

       cout << "---------\nt = " << t << endl;
       cout << "tau = " << tau << endl;

       k1 = solver.assembleRHS(lhs);

       solver.alphaNext = solver.alphaPrev + k1 * tau;
       lhs = solver.correctNonOrtho(solver.alphaNext);

       limiter.limit(lhs);


       k2 = solver.assembleRHS(lhs);

       solver.alphaNext = solver.alphaPrev + (k1 + k2) * 0.5 * tau;
       lhs = solver.correctNonOrtho(solver.alphaNext);

       //cout << "before limiting" << solver.alphaNext[49] << endl;

       limiter.limit(lhs);

       //cout << "after limiting" << solver.alphaNext[49] << endl;

       if (iT % freqWrite == 0)
       {
           //string fileName = "alphaCoeffs/" + to_string((long double)t);
           string fileName = "alphaCoeffs/" + to_string(t);

           ofstream output;
           output.open(fileName + ".vtk");

           mesh.exportMeshVTK(output);

           output << "CELL_DATA " << mesh.nCells << endl;

           output << "FIELD rho 1" << endl;
           output << "cellscalar" << " 1" << " " << mesh.nCells << " float" << endl;

           for (int i = 0; i < mesh.nCells; ++i)
              output << mesh.cells[i]->reconstructSolution(mesh.cells[i]->getCellCenter(),0) << endl;

           output << "FIELD rhoU 1" << endl;
           output << "rhoU " << " 3" << " " << mesh.nCells << " float" << endl;

           for (int i = 0; i < mesh.nCells; ++i)
           {
              output << mesh.cells[i]->reconstructSolution(mesh.cells[i]->getCellCenter(),1) << ' ';
              output << mesh.cells[i]->reconstructSolution(mesh.cells[i]->getCellCenter(),2) << ' ';
              output << 0.0 << endl;
           }

           output << "FIELD e 1" << endl;
           output << "e" << " 1" << " " << mesh.nCells << " float" << endl;

           for (int i = 0; i < mesh.nCells; ++i)
              output << mesh.cells[i]->reconstructSolution(mesh.cells[i]->getCellCenter(),4) << endl;

           output << "POINT_DATA " << mesh.nodes.size() << endl;

           output.close();
       }

       // get limited "lhs"
       solver.alphaNext = solver.correctPrevIter(lhs);
       solver.alphaPrev = solver.alphaNext;

       iT++;

       t2 = clock();

       dynamicTimeController.updateTimeStep();

       tau = dynamicTimeController.getNewTau();

       cout << "step time: " << (float)(t2 - t1) / CLOCKS_PER_SEC << endl;
    }

    cout << "=========\nElapsed time = " << (float)(t2 - t00) / CLOCKS_PER_SEC << endl;
    cout << "---------\nEND \n";

    //cin.get();
	return 0;
}

