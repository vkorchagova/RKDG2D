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
#include "IndicatorKXRCF.h"

using namespace std;


int main(int argc, char** argv)
{    
    // Mesh parameters

    double Lx = 4.0;
    double Ly = 4.0;

    int nx = 20;
    int ny = 40;

    // foolish left boundary!!!!

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


    //cout << "Indicator test..." << endl;
    //vector<double> ind;

	//Testing Indicator
    //IndicatorKXRCF indicator(mesh);
    //auto ind = indicator.checkDiscontinuities();
    
    //for (size_t q = 0; q < ind.size(); ++q)
	//cout << "q = " << q << ", ind = " << ind[q] << endl;


    // time cycle paramentes

    double Co = 0.1;
    double tEnd = 2.01;

    double tau = min(mesh.cells[0]->h().x(),mesh.cells[0]->h().y()) * Co;
        // sound speed = 1 --- const in acoustic problems
        // only for uniform mesh hx and hy are similar for all cells


    // run Runge --- Kutta 2 TVD

    vector<numvector<double, 5*nShapes>> k1, k2;

    k1.resize(mesh.nCells);
    k2.resize(mesh.nCells);

    clock_t t1, t2, t00;

    t00 = clock();

    for (double t = tau; t < tEnd; t += tau)
    {
        t1 = clock();
       //string fileName = "alphaCoeffs/" + to_string((long double)t);
       string fileName = "alphaCoeffs/" + to_string(t);

       ofstream output;
       output.open(fileName);

       cout << "---------\nt = " << t << endl;

       k1 = solver.assembleRHS(solver.alphaPrev);
       solver.alphaNext = solver.alphaPrev + k1 * tau;

       //problem.setAlpha(solver.alphaNext);

       /*
	   ind = indicator.checkDiscontinuities();

	   for (size_t i = 0; i < mesh.nCells; ++i)
	   {
		   if (ind[i] > 1.0)
		   {
			   for (int j = 0; j < 5; ++j)
			   {
				   solver.alphaNext[i][j*nShapes + 1] = 0.0;
				   solver.alphaNext[i][j*nShapes + 2] = 0.0;
			   }
		   }
	   }
*/

       k2 = solver.assembleRHS(solver.alphaNext);
       solver.alphaNext = solver.alphaPrev + (k1 + k2) * 0.5 * tau;

       //problem.setAlpha(solver.alphaNext);

       // limit solution


       solver.write(output,solver.alphaNext);

       solver.alphaPrev = solver.alphaNext;

       output.close();
       t2 = clock();

       cout << "step time: " << (float)(t2 - t1) / CLOCKS_PER_SEC << endl;
    }

    cout << "---------\nElapsed time = " << (float)(t2 - t00) / CLOCKS_PER_SEC << endl;

    cout << "END \n";

//    int aaa;
//    cin >> aaa;

	cin.get();
	return 0;
}

