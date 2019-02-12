#pragma once

#include <mpi.h>
#include <omp.h>

//- Number of conservative variables
static const int PhysDim = 5;
static const int dimPh = 5;

//- Number of basis functions
static const int nShapes = 3;

// Coeffs vector size
static const int dimS = dimPh * nShapes;

const int NumThreads = 4; // Number of threads for teh OpenMP

//- Initialisaton of tools for computations
//void Initialize() {}

//- list of all variables 
enum Variables 
{ 
    RHO  = 0,  
    RHOU = 1,  
    RHOV = 2, 
    RHOW = 3, 
    E    = 4
}; 

//- list if initial cases 
enum CaseInit 
{ 
    SodX, SodY, SodDiag, SodCircle, Blast, Const, Acoustic, Acoustic1D, ForwardStep, Sedov
}; 

//- list of BC types
enum CaseBound 
{ 
    Inf, Free, Diag, Wall, PeriodicSquare 
};

//- List of global MPI variables



