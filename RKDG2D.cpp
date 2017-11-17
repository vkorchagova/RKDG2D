// RKDG 2D v.0.1
// Structured rectangular mesh


#include <stdio.h>
#include <iostream>
#include "Mesh2D.h"

using namespace std;



int main(int argc, char** argv)
{
    int nx = 4;
    int ny = 5;
    double Lx = 4;
    double Ly = 4;

    Mesh2D mesh(nx, ny, Lx, Ly);

	//for (int i = 0; i < mesh.neighbCellsVerEdges.size(); ++i)
	//	cout << mesh.neighbCellsVerEdges[i] << endl;

	cout << mesh.globalToLocal(0, {1,0.8});
	cout << mesh.localToGlobal(0, {1,1});

	int hz;
	cin >> hz;


	return 0;
}

