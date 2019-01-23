// DG_MPI.cpp
//


#if !defined(__linux__)
#include <direct.h>
#endif

#include <sys/stat.h>


#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>
#include <time.h>
#include "numvector.h"

#include "defs.h"		//- Basic arithmetics
#include "Mesh.h"		//- The mesh
#include "Basis.h"		//- All about the basis functions for DG
#include "compService.h"//- Integration
#include "Physics.h"	//- The physical models
#include "Solution.h"   //- Solution storage
#include "Writer.h"

#include <omp.h>

using namespace std;


void setDefinedCoefficients(string fileName, const Mesh& M, Solution& sln)
{
    int nCells = M.cells.size();
    sln.SOL.resize(nCells);

    ifstream reader;
    reader.open(fileName);

    if (!reader.is_open())
    {
        cout << "File " << fileName << " is not found\n";
        exit(0);
    }// if open

    numvector<double, dimS> rhs;

    for (int k = 0; k < nCells; ++k)
    {
        for (int j = 0; j < dimS; ++j)
            reader >> sln.SOL[k][j];
    }// for k

    reader.close();
}


int main(int argc, char* argv[])
{

    //-----------------------

    string meshFileName = "mesh2D";
    string solFileName = "";

    cout << "Before mesh..." << endl;

    Mesh mesh(meshFileName);

    cout << "Mesh OK" << endl;
    Basis basis(mesh.cells);
    Solution solution(basis);    
    Physics physics;
    physics.cpcv = 1.4;     // TODO-TODO-TODO...

    Writer writer(mesh, solution, physics);
    
    ifstream timePointsFile("times");

    if (!timePointsFile.is_open())
    {
        cout << "File " << "times" << " is not found\n";
        exit(0);
    };

    while (timePointsFile.peek() != EOF)
    {
        timePointsFile >> solFileName;

        cout << "Proceeding frame " << solFileName << "..." << endl;

        setDefinedCoefficients("alphaCoeffs/" + solFileName + ".dat", mesh, solution);
 
        writer.exportFrameVTK("alphaCoeffs/" + solFileName + ".vtk");
    }

    timePointsFile.close();
    

    cout << "THE END" << endl;
	
	return 0;
}
