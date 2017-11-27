//- RKDG 2D v.0.1
//  Structured rectangular mesh


#include <stdio.h>
#include <iostream>
#include "Mesh2D.h"
#include "Problem.h"

using namespace std;



int main(int argc, char** argv)
{
    int nx = 4;
    int ny = 5;
    double Lx = 4;
    double Ly = 4;

    // Get mesh
    Mesh2D mesh(nx, ny, Lx, Ly);

    // Get solver
    Problem problem(mesh);

	problem.setInitialConditions();

	for (int i = 0; i < mesh.cells.size(); ++i)
		cout << problem.alphaPrev[i] << endl;

	int aaa;
	cin >> aaa;


	return 0;
}

