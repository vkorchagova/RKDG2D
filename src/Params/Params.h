#pragma once

#include <mpi.h>
#include <omp.h>

/// Number of conservative variables
static const int dimPh = 5; //

/// Number of basis functions
static const int nShapes = 3;

// Number of mean fields (mean Ux + mean Uy + mean p)
static const int dimMean = 3;

/// Number of exported variables
static const int dimExp = 6 + dimMean; // 5 conservative + pressure + mean Ux + mean Uy + mean p

// Coeffs vector size
static const int dimS = dimPh * nShapes;


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
    SodX, SodXCovol, SodY, SodDiag, SodCircle, Blast, Const, Acoustic, Acoustic1D, ForwardStep, Sedov, DoubleMach, Ladenburg, LadenburgRibbon, ShuOsher, AstroTest, BigPulse, acousticPulse, Munday, ConvDivNozzle
}; 

/// list of BC types
enum CaseBound 
{ 
    Inf, Free, Diag, Wall, PeriodicSquare 
};



