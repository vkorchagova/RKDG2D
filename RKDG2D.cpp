//- RKDG 2D v.0.1
//  Structured rectangular mesh



#include <stdio.h>
#include <iostream>
#include <string>
#include <time.h>
#include "defs.h"
#include "Mesh2D.h"
#include "Solver.h"
#include "FluxLLF.h"
#include "FluxHLL.h"
#include "FluxHLLC.h"
#include "IndicatorEverywhere.h"
#include "IndicatorKXRCF.h"
#include "LimiterFinDiff.h"
#include "LimiterMUSCL.h"
#include "LimiterWENOS.h"

using namespace std;


int main(int argc, char** argv)
{    
    // Mesh parameters

    double Lx = 1.0;
    double Ly = 1.0;

    int nx = 40;
    int ny = 40;

//    double Lx = 4.0;
//    double Ly = 4.0;

//    int nx = 20;
//    int ny = 20;

    // Time parameters

    double Co = 0.1;
    double tEnd = 0.2;

    // ---------------

    // Initialize mesh
    Mesh2D mesh(nx, ny, Lx, Ly);

    mesh.exportMesh();

    // Initialize problem
    Problem problem;

    // Initialize flux
    FluxLLF numFlux(problem);

    // Initialize solver
    Solver solver(mesh, problem, numFlux);

    // Initialize indicator
    IndicatorKXRCF indicator(mesh);

    //Initialize limiter
    LimiterWENOS limiter(indicator,problem);

    // ---------------

    // Set initial conditions
    solver.setInitialConditions();

    // Set boundary conditions
    solver.setBoundaryConditions();

    // time step

    double tau = min(mesh.cells[0]->h().x(),mesh.cells[0]->h().y()) * Co;
        // sound speed = 1 --- const in acoustic problems
        // only for uniform mesh hx and hy are similar for all cells


    // run Runge --- Kutta 2 TVD

    vector<numvector<double, 5*nShapes>> k1, k2;

    k1.resize(mesh.nCells);
    k2.resize(mesh.nCells);

    clock_t t1, t2, t00;

    t00 = clock();

    for (double t = tau; t <= tEnd + 0.5*tau; t += tau)
    {
        t1 = clock();
       //string fileName = "alphaCoeffs/" + to_string((long double)t);
       string fileName = "alphaCoeffs/" + to_string(t);

       ofstream output;
       output.open(fileName);

       cout << "---------\nt = " << t << endl;

       k1 = solver.assembleRHS(solver.alphaPrev);
       solver.alphaNext = solver.alphaPrev + k1 * tau;

       limiter.limit(solver.alphaNext);

       k2 = solver.assembleRHS(solver.alphaNext);
       solver.alphaNext = solver.alphaPrev + (k1 + k2) * 0.5 * tau;

       limiter.limit(solver.alphaNext);

       solver.write(output,solver.alphaNext);

       solver.alphaPrev = solver.alphaNext;

       output.close();
       t2 = clock();

       cout << "step time: " << (float)(t2 - t1) / CLOCKS_PER_SEC << endl;
    }

    cout << "=========\nElapsed time = " << (float)(t2 - t00) / CLOCKS_PER_SEC << endl;
    cout << "---------\nEND \n";

    //cin.get();
	return 0;
}

