//- RKDG 2D v.0.1
//  Structured rectangular mesh


#include <stdio.h>
#include <iostream>

#include <sstream>
#include "Mesh2D.h"
#include "Problem.h"
#include "fluxllf.h"
#include "defs.h"


int main(int argc, char** argv)
{
    // Mesh parameters
    int nx = 20;
    int ny = 20;
    double Lx = 4.0;
    double Ly = 4.0;

    // Courant number
    double Co = 0.25;

    // Get mesh
    Mesh2D mesh(nx, ny, Lx, Ly);

    mesh.exportMesh();

    // Get solver
    Problem problem(mesh);

    problem.setInitialConditions();
    problem.applyBoundary(problem.alphaPrev);


    //Get flux
    FluxLLF numFlux (problem);


    // time cycle

    // TODO: dynamic time step!!!



    double tau = min(mesh.hx,mesh.hy) * Co; // sound speed = 1 --- const in acoustic problems

    double tEnd = 2.01;

    vector<numvector<double, 5*nShapes>> k1, k2;

    int nCells = mesh.nInternalCells + mesh.nGhostCells ;

    k1.resize(nCells);
    k2.resize(nCells);

    // open ofstream for coeffs
    string fileName = "alphaCoeffs/" + to_string(0.0);
    ofstream output;
    output.open(fileName);

    problem.write(output,problem.alphaPrev);

    output.close();

    // run Runge --- Kutta 2 TVD
    for (double t = tau; t < tEnd; t += tau)
    {
        string fileName = "alphaCoeffs/" + to_string(t);

        ofstream output;
        output.open(fileName);

        cout << "t = " << t << endl;

        k1 = numFlux.getRHS(problem.alphaPrev);
        problem.alphaNext = problem.alphaPrev + k1 * tau;

        problem.applyBoundary(problem.alphaNext);

        k2 = numFlux.getRHS(problem.alphaNext);
        problem.alphaNext = problem.alphaPrev + (k1 + k2) * 0.5 * tau;

        problem.applyBoundary(problem.alphaNext);

        problem.write(output,problem.alphaNext);

        problem.alphaPrev = problem.alphaNext;

        output.close();
    }

    cout << "END \n";

    //int aaa;
    //cin >> aaa;


	return 0;
}

