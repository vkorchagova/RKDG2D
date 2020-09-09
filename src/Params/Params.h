#pragma once

#include <mpi.h>
#include <omp.h>

/// Number of conservative variables
static const int dimPh = 5; //

/// Number of basis functions
static const int nShapes = 3;

/// Number of exported variables
static const int dimExp = 6; // conservative + pressure

// Coeffs vector size
static const int dimS = dimPh * nShapes;

/// Number of grad-components
static const int dimGrad = 8;
static const int dimGradFull = 15;

/// Coeffs vector size for grad
static const int dimGradCoeff = dimGrad * nShapes;

/// Initialisaton of tools for computations
//void Initialize() {}

/// list of all variables 
enum Variables 
{ 
    RHO  = 0,  
    RHOU = 1,  // and U
    RHOV = 2,  // and V
    RHOW = 3,  // and W
    E    = 4,
    P    = 5
}; 

/// list if initial cases 
enum CaseInit 
{ 
    SodX, SodXCovol, SodY, SodDiag, SodCircle, Blast, Const, Acoustic, Acoustic1D, ForwardStep, Sedov, DoubleMach, Ladenburg, LadenburgRibbon, ShuOsher, AstroTest, Blasius, Vortex, CylinderFlow
}; 

/// list of BC types
enum CaseBound 
{ 
    Inf, Free, Diag, Wall, PeriodicSquare 
};



