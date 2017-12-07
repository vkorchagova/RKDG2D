//- RKDG 2D v.0.1
//  Structured rectangular mesh


#include <stdio.h>
#include <iostream>
#include <string>
#include <sstream>
#include "Mesh2D.h"
#include "Problem.h"
#include "fluxllf.h"
#include "defs.h"

using namespace std;



int main(int argc, char** argv)
{
    // Mesh parameters
    int nx = 25;
    int ny = 20;
    double Lx = 5;
    double Ly = 4;

    // Courant number
    double Co = 0.25;

    // Get mesh
    Mesh2D mesh(nx, ny, Lx, Ly);

    mesh.exportMesh();

    // Get solver
    Problem problem(mesh);

    problem.setInitialConditions();
    problem.applyBoundary(problem.alphaPrev);

    //for (int i = 0; i < mesh.cells.size(); ++i)
      //  cout << i << ' ' << problem.alphaPrev[i] << endl;

    //Get flux
    FluxLLF numFlux (problem);

    //numvector<double,2> pt = {1.0,0.4};
    //numvector<double,5> res,s1,s2;
    //res=problem.reconstructSolution(problem.alphaPrev[0],pt,0);
    //s1=problem.reconstructSolution(problem.alphaPrev[0],pt,0);
    //s2=problem.reconstructSolution(problem.alphaPrev[1],pt,1);
    //res=problem.fluxG(problem.reconstructSolution(problem.alphaPrev[0],pt,0));
    //for(int p=0;p<5;++p)
     //   cout<<res[p]<<endl;
    //cout<< problem.c_av(s1,s2) << endl;



    //cout << numFlux.problem->mesh->cellCenters.size() << endl;
    //problem.write(cout,problem.alphaPrev[9]);

    //vector<numvector<double, 5*nShapes>> rhs = numFlux.getRHS(problem.alphaPrev);

    // time cycle

    // TODO: dynamic time step!!!



    double tau = min(mesh.hx,mesh.hy) * Co; // sound speed = 1 --- const in acoustic problems


    double tEnd = 2.01;

    vector<numvector<double, 5*nShapes>> k1, k2;

    int nCells = mesh.nInternalCells + mesh.nGhostCells ;

    k1.resize(nCells);
    k2.resize(nCells);

    string fileName = "alphaCoeffs/" + to_string(0.0);
    ofstream output;
    output.open(fileName);

    problem.write(output,problem.alphaPrev);

    output.close();


/*
    numvector<double, 2> point = {mesh.hx * 0.5, 0.0};
    numvector<double,5> temp;
    temp = problem.reconstructSolution(problem.alphaNext[nx*ny],point,nx*ny);
    cout << "comment 01_" << temp[0] << endl;

    //if (i==0 &&j==0) cout << "comment " << temp[0] << endl;
*/

/*
    //problem.write(cout,problem.alphaPrev);
    k1 = numFlux.getRHS(problem.alphaPrev);
    //problem.write(cout,k1);
    problem.alphaNext = problem.alphaPrev + k1 * tau;
    problem.applyBoundary(problem.alphaNext);
    //problem.write(cout,problem.alphaNext);
    k2 = numFlux.getRHS(problem.alphaNext);
    //problem.write(cout,k2);
    problem.alphaNext = problem.alphaPrev + (k1 + k2) * 0.5 * tau;
    problem.applyBoundary(problem.alphaNext);
    problem.alphaPrev = problem.alphaNext;
    //problem.write(cout,problem.alphaPrev);


    k1 = numFlux.getRHS(problem.alphaPrev);
    //problem.write(cout,k1);
    problem.alphaNext = problem.alphaPrev + k1 * tau;
    problem.applyBoundary(problem.alphaNext);
    //problem.write(cout,problem.alphaNext);
    k2 = numFlux.getRHS(problem.alphaNext);
    //problem.write(cout,k2);
    problem.alphaNext = problem.alphaPrev + (k1 + k2) * 0.5 * tau;
    problem.applyBoundary(problem.alphaNext);
    problem.alphaPrev = problem.alphaNext;
    problem.write(cout,problem.alphaPrev);

*/

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

    problem.write(cout,problem.alphaPrev);


    cout << "END \n";

    //int aaa;
    //cin >> aaa;


	return 0;
}

