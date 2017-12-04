#include "Flux.h"

using namespace std;

// ------------------ Constructors & Destructors ----------------


Flux::Flux(Problem &prb)
{
    problem = &prb;
    mesh = problem->mesh;
}


Flux::~Flux()
{
}

// ------------------ Private class methods --------------------

// ------------------ Public class methods ---------------------

