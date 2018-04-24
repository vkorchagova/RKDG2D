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
    double tEnd = 0.2e-6;

    double initDeltaT = 5e-3;
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
    FluxHLL numFlux(problem);

    // Initialize solver
    Solver solver(mesh, problem, numFlux);

    // Initialize time controller
    TimeControl dynamicTimeController(mesh,Co,maxDeltaT,maxTauGrowth,initDeltaT,isDynamicTimeStep);

    // Initialize indicator
    IndicatorKXRCF indicator(mesh, problem);

    // Initialize limiter
    LimiterFinDiff limiter(indicator, problem);

    // ---------------

    #if !defined(__linux__)
        _mkdir("alphaCoeffs");
    #else
        mkdir("alphaCoeffs", S_IRWXU | S_IRGRP | S_IROTH);
    #endif

    // Set initial conditions
    solver.setInitialConditions();

    limiter.limit(solver.alphaPrev);

    ofstream output;
    output.open("alphaCoeffs/0.000000");

    solver.write(output,solver.alphaPrev);

    output.close();

//    for (const shared_ptr<Cell> cell : mesh.cells)
//    {
//        cout << "neibs for cell #" << cell-> number << ": ";
//        for (const shared_ptr<Cell> neib : cell->neibCells)
//            cout << neib->number << ' ';
//        cout << endl;
//    }

   // time step

    double tau = initDeltaT;


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

       k1 = solver.assembleRHS(solver.alphaPrev);

       solver.alphaNext = solver.alphaPrev + k1 * tau;
       solver.correctNonOrtho(solver.alphaNext);

       limiter.limit(solver.alphaNext);


       k2 = solver.assembleRHS(solver.alphaNext);

       solver.alphaNext = solver.alphaPrev + (k1 + k2) * 0.5 * tau;
       solver.correctNonOrtho(solver.alphaNext);

       //cout << "before limiting" << solver.alphaNext[49] << endl;

       limiter.limit(solver.alphaNext);

       //cout << "after limiting" << solver.alphaNext[49] << endl;

       if (iT % freqWrite == 0)
       {
           //string fileName = "alphaCoeffs/" + to_string((long double)t);
           string fileName = "alphaCoeffs/" + to_string(t);

           ofstream output;
           output.open(fileName);

           solver.write(output,solver.alphaNext);

           output.close();
       }

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

