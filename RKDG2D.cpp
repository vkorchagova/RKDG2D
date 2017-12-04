//- RKDG 2D v.0.1
//  Structured rectangular mesh


#include <stdio.h>
#include <iostream>
#include "Mesh2D.h"
//#include "Problem.h"
//#include "fluxllf.h"

using namespace std;



int main(int argc, char** argv)
{
    int nx = 40;
    int ny = 50;
    double Lx = 4;
    double Ly = 4;

    // Get mesh
    Mesh2D mesh(nx, ny, Lx, Ly);

  //  mesh.exportMesh();

    // Get solver
    //Problem problem(mesh);



    //problem.setInitialConditions();

    //for (int i = 0; i < mesh.cells.size(); ++i)
      //  cout << i << ' ' << problem.alphaPrev[i] << endl;

    //Get flux
    //FluxLLF numFlux (problem);

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

    //problem.write(cout,rhs[9]);

    cout << "END \n";

    int aaa;
    cin >> aaa;


	return 0;
}

