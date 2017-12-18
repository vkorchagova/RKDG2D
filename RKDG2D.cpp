//- RKDG 2D v.0.1
//  Structured rectangular mesh



#include <stdio.h>
#include <iostream>
#include "defs.h"
#include "Mesh2D.h"
#include "Solver.h"
#include "FluxLLF.h"

using namespace std;




int main(int argc, char** argv)
{    
    // Mesh parameters

    double Lx = 8.0;
    double Ly = 8.0;

    int nx = 20;
    int ny = 40;

    // Initialize mesh
    Mesh2D mesh(nx, ny, Lx, Ly);

    mesh.exportMesh();

    // Initialize problem
    Problem problem;

    // Initialize flux
    FluxLLF numFlux (problem);

    // Initialize solver
    Solver solver(mesh, problem, numFlux);

    // Set initial conditions
    solver.setInitialConditions();

    // Set boundary conditions
    solver.initBoundaryConditions();


//    cout << "Init cond:\n " ;
//    for (int i = 0; i < mesh.cells.size(); ++i)
//        cout << problem.alpha[i] << endl;

    // time cycle paramentes

    double Co = 0.25;
    double tEnd = 2.01;

    double tau = min(mesh.cells[0]->h().x(),mesh.cells[0]->h().y()) * Co;
        // sound speed = 1 --- const in acoustic problems
        // only for uniform mesh hx and hy are similar for all cells


    // run Runge --- Kutta 2 TVD

    vector<numvector<double, 5*nShapes>> k1, k2;

    k1.resize(mesh.nCells);
    k2.resize(mesh.nCells);

    for (double t = tau; t < tEnd; t += tau)
    {
       string fileName = "alphaCoeffs/" + to_string(t);

       ofstream output;
       output.open(fileName);

       cout << "t = " << t << endl;

       k1 = solver.assembleRHS(solver.alphaPrev);
       solver.alphaNext = solver.alphaPrev + k1 * tau;

       k2 = solver.assembleRHS(solver.alphaNext);
       solver.alphaNext = solver.alphaPrev + (k1 + k2) * 0.5 * tau;

       solver.write(output,solver.alphaNext);

       solver.alphaPrev = solver.alphaNext;

       output.close();
    }

    cout << "END \n";

//    int aaa;
//    cin >> aaa;


	return 0;
}

