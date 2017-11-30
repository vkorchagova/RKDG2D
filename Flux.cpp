#include "Flux.h"

using namespace std;

// ------------------ Constructors & Destructors ----------------


Flux::Flux(Problem &prb)
{
    problem = &prb;
}


Flux::~Flux()
{
}

// ------------------ Private class methods --------------------

// ------------------ Public class methods ---------------------

vector<numvector<double,5*nShapes>> Flux::getRHS()
{
    double nEdgesHor = problem->mesh->edgesHor.size();
    double nEdgesVer = problem->mesh->edgesVer.size();

    vector<double> integrateFluxesHor(nEdgesHor);
    vector<double> integrateFluxesVer(nEdgesVer);

    GaussIntegrator GP;

    for (int i = 0; i < nEdgesHor; ++i)
    {
        numvector<double,5> res;

        for (int j = 0; j < 5; ++j)
        {
            double a = problem->mesh->nodes[problem->mesh->edgesHor[i][0]][0];
            double b = problem->mesh->nodes[problem->mesh->edgesHor[i][1]][0];
            double y = problem->mesh->nodes[problem->mesh->edgesHor[i][0]][1];

            //function<double(double)> f = [](double(double))
           // {
             //   return 0;
            //}
        }
    }


}
