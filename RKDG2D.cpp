//- RKDG 2D v.0.2
//  Unstructured 2D mesh

#if !defined(__linux__)
#include <direct.h>
#endif

#include <sys/stat.h>

#include <stdio.h>
#include <iostream>
#include <string>
#include <omp.h>

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
#include "LimiterRiemannWENOS.h"


using namespace std;



int main(int argc, char** argv)
{    
    // Problem

    string caseName = "monopole";

    omp_set_num_threads(1);

    // Time parameters

    double tStart = 0.0;
    double tEnd = 0.005;

    double outputInterval = 0.005;
    double initDeltaT = 5e-3;

    bool   defCoeffs = false; // true if alpha coefficients for start time are defined

    bool   isDynamicTimeStep = false;
    double Co = 0.5;
    double maxDeltaT = 1.0;
    double maxTauGrowth = 1.1;

    // ---------------

    // Initialize time
    Time time;

    // Initialize problem
    Problem problem(caseName,time);

    // Initialize mesh
    Mesh2D mesh("mesh2D",problem);

    // Set BC
    problem.setBoundaryConditions(caseName, mesh.patches);

    // Initialize flux
    FluxHLLC numFlux(problem);

    // Initialize solver
    Solver solver(mesh, problem, numFlux);

    // Initialize time controller
    TimeControl dynamicTimeController(mesh,Co,maxDeltaT,maxTauGrowth,initDeltaT,isDynamicTimeStep);

    // Initialize indicator
    IndicatorNowhere indicator(mesh, problem);

    // Initialize limiter
    LimiterRiemannWENOS limiter(indicator, problem);

    // ---------------

    #if !defined(__linux__)
        _mkdir("alphaCoeffs");

    #else
        mkdir("alphaCoeffs", S_IRWXU | S_IRGRP | S_IROTH);
    #endif

    cout << "---------\nt = " << tStart << endl;

    // Set initial conditions
    if (defCoeffs)
    {
        solver.setDefinedCoefficients("alphaCoeffs/" + to_string(tStart));
    }
    else
    {
        solver.setInitialConditions();
        limiter.limit(solver.alphaPrev);
        solver.writeSolutionVTK("alphaCoeffs/sol_" + to_string(tStart));
    }


    // time step

    double tau = initDeltaT;

    vector<numvector<double, 5*nShapes>> lhs = solver.alphaPrev; // coeffs
    //vector<numvector<double, 5*nShapes>> lhsPrev = lhs;

    // run Runge --- Kutta 3 TVD

    vector<numvector<double, 5*nShapes>> k1, k2, k3;

    k1.resize(mesh.nCells);
    k2.resize(mesh.nCells);
    k3.resize(mesh.nCells);

    double t1, t2, t00;

    t00 = omp_get_wtime();

    int iT = 1; //iteration number
    int nOutputSteps = 1;

    for (double t = tStart + tau; t <= tEnd + 0.5*tau; t += tau)
    {
       t1 = omp_get_wtime();

       time.updateTime(t);

       cout.precision(6);
       cout << "---------\nt = " << t << endl;
       cout << "tau = " << tau << endl;

       k1 = solver.assembleRHS(lhs);

       solver.alphaNext = solver.alphaPrev + k1 * tau;
       lhs = solver.correctNonOrtho(solver.alphaNext);

       limiter.limit(lhs);
       solver.alphaNext = solver.correctPrevIter(lhs);

//       solver.writeSolutionVTK("alphaCoeffs/sol_" + to_string(t)+"RK1");
       //solver.write("alphaCoeffs/" + to_string(t)+"RK1",lhs);

//       double totalEnergy = 0.0;

//       for (const shared_ptr<Cell> cell : mesh.cells)
//       {

//           function<double(const Point&)> eTotal = [&](const Point& x)
//           {
//               return cell->reconstructSolution(x,4);
//           };

//           totalEnergy += cell->integrate(eTotal);
//       }
//       cout.precision(16);
//       cout << "total energy = " << totalEnergy << endl;


       k2 = solver.assembleRHS(lhs);

       //solver.alphaNext = solver.alphaPrev + (k1 + k2) * 0.5 * tau;
       solver.alphaNext = solver.alphaPrev*0.75 + solver.alphaNext*0.25 + k2*0.25*tau;
       lhs = solver.correctNonOrtho(solver.alphaNext);

       limiter.limit(lhs);
       solver.alphaNext = solver.correctPrevIter(lhs);

//       totalEnergy = 0.0;

//       for (const shared_ptr<Cell> cell : mesh.cells)
//       {

//           function<double(const Point&)> eTotal = [&](const Point& x)
//           {
//               return cell->reconstructSolution(x,4);
//           };

//           totalEnergy += cell->integrate(eTotal);
//       }
//       cout.precision(16);
//       cout << "total energy = " << totalEnergy << endl;

       k3 = solver.assembleRHS(lhs);
       double i13 = 0.3333333333333333;

       solver.alphaNext = solver.alphaPrev*i13 + solver.alphaNext*2*i13 + k3*2*i13*tau;
       lhs = solver.correctNonOrtho(solver.alphaNext);

       limiter.limit(lhs);
       solver.alphaNext = solver.correctPrevIter(lhs);


       // check energy conservation
       double totalEnergy = 0.0;

#pragma omp parallel for 
       //for (const shared_ptr<Cell> cell : mesh.cells)
       for ( size_t i = 0; i < mesh.cells.size(); ++i )
       {
          const shared_ptr<Cell> cell = mesh.cells[i];
           function<double(const Point&)> eTotal = [&](const Point& x)
           {
               return cell->reconstructSolution(x,4);
           };

           totalEnergy += cell->integrate(eTotal);
       }
       cout.precision(16);
       cout << "total energy = " << totalEnergy << endl;

       //cout << "after limiting" << solver.alphaNext[49] << endl;
//       solver.write("alphaCoeffs/" + to_string(t)+"RK2bl",lhs);

       //if (iT % freqWrite == 0)
       // if (fabs(t - nOutputSteps * outputInterval) < 1e-10)
       // {
       //     //string fileName = "alphaCoeffs/" + to_string((long double)t);

       //     solver.writeSolutionVTK("alphaCoeffs/sol_" + to_string(t));
       //     solver.write("alphaCoeffs/" + to_string(t),lhs);
       //     nOutputSteps++;
       // }



       // check local internal energy balance

//       double internalEnergyB = 0.0;

//       {
//           const shared_ptr<Cell> cell = mesh.cells[50];

//           function<double(const Point&)> eIn = [&](const Point& x)
//           {
//               numvector<double, 5> solPrev  = cell->reconstructSolution(x, lhsPrev[50]);
//               numvector<double, 5> solNext  = cell->reconstructSolution(x, lhs[50]);

//               double eInPrev = solPrev[4] - 0.5 * (sqr(solPrev[1]) + sqr(solPrev[2]) / solPrev[0]);
//               double eInNext = solNext[4] - 0.5 * (sqr(solNext[1]) + sqr(solNext[2]) / solNext[0]);

//               double ux = lhsPrev[cell->number][1 * nShapes + 1] * cell->offsetPhi[1] / solPrev[0];
//               double vy = lhsPrev[cell->number][2 * nShapes + 2] * cell->offsetPhi[2] / solPrev[0];
//               double pDivU = problem.getPressure(solPrev) * (ux + vy);

//               return eInNext - eInPrev - pDivU;
//           };

//          internalEnergyB = cell->integrate(eIn);
//       }




       //cout << "internal energy balance = " << internalEnergyB << endl;




       solver.alphaPrev = solver.alphaNext;
       //lhsPrev = lhs;

       //iT++; 

       dynamicTimeController.updateTimeStep();

       //tau = min(nOutputSteps * outputInterval - t, dynamicTimeController.getNewTau());
       tau = dynamicTimeController.getNewTau();

       t2 = omp_get_wtime();

       cout << "step time: " << t2 - t1 << endl;
    }

    cout << "=========\nElapsed time = " << t2 - t00 << endl;
    cout << "---------\nEND \n";

    //cin.get();
	return 0;
}

