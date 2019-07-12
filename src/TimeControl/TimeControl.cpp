#include "TimeControl.h"
#include <algorithm>
#include <iomanip>

using namespace std;



TimeControl::TimeControl(
        const Mesh& msh,
        //const Problem& prb, 
        const double tStart, 
        const double tEnd, 
        const double initTau, 
        const double outputInterval,
        const string listingPath
    ) : 
    M(msh), 
    //problem(prb),
    t(tStart), 
    tau(initTau), 
    tEnd(tEnd), 
    outputInterval(outputInterval)
{
    tauNew = tau;
    tOld = t;
    outputTime = tStart + outputInterval;
    timeListing.open(listingPath.c_str());

    if (!timeListing.is_open())
    {
        cout << "File " << listingPath << " could not be opened\n";
        exit(0);
    }

    timeListing << setprecision(6);
    timeListing << 0.000000 << endl;
    //write the construction with initialization from Params.h or some other data source;
}

TimeControl::~TimeControl()
{
    timeListing.close();
}

void TimeControl::updateTimeStep(double MSpeed) // maxSpeed should be defined on each cell
{
    if (isDynamic) // template for the code of Vik
    {
        /*
        double factCo = 0.0;
        double relTau = 0.0;

        vector<double> newTauLocal;
        newTauLocal.reserve(M.nRealCells);

        vector<double> maxUL;
        maxUL.reserve(M.nRealCells);

//#pragma omp parallel for \
//shared (mesh) \
//default (none)
        // get maxUL for bound edges 
        for (const shared_ptr<Boundary>& bcond : prb.bc)
        {
            //#pragma omp parallel for \
                shared(myRank, nGP, numFluxes, bcond) \
                private (solLeft, solRight, gpFluxes) \
                default(none)
            //for (const shared_ptr<Edge>& e : bcond->patch.edgeGroup)
            for (int iEdge = 0; iEdge < bcond->patch.edgeGroup.size(); ++iEdge)
            {
                const shared_ptr<Edge>& e = bcond->patch.edgeGroup[iEdge];
                int iCellLeft  = e->neibCells[0]->number;
                double uMaxOneEdge = 0.0;


                for (size_t i = 0; i < 2; ++i)
                {
                    Point& node = e->nodes[i];
                    Point& eNormal = e->n;

                    solLeft  = rotate(sln.reconstruct(iCellLeft,  node), eNormal);
                    solRight = bcond->getSolOuter(solLeft);

                    numvector<double, 5> lambda = problem.lambdaF(solLeft,solRight);
                    uMaxOneEdge = max( max(fabs(lambda[0]), fabs(lambda[4])), uMaxOneEdge);
                }// for GP

                maxUL[e->number] = uMaxOneEdge * e->length();
            }// for bound edges
        } // for bconds 

        // get maxUL for internal edges
        for (int iEdge = M.nEdgesBound; iEdge < M.nRealEdges; ++iEdge)
        {
            const shared_ptr<Edge>& e = M.edges[iEdge];

            double uMaxOneEdge = 0.0;

            int iCellLeft  = e->neibCells[0]->number;
            int iCellRight = e->neibCells[1]->number;

            for (int i = 0; i < 2; ++i)
            {
                const Point& node = e->nodes[i];
                const Point& eNormal = e->n;
                
                solLeft  = rotate(sln.reconstruct(iCellLeft,  node), eNormal);
                solRight = rotate(sln.reconstruct(iCellRight, node), eNormal);

                numvector<double, 5> lambda = problem.lambdaF(solLeft,solRight);
                uMaxOneEdge = max( max(fabs(lambda[0]), fabs(lambda[4])), uMaxOneEdge);

                gpFluxes[iGP] = inverseRotate(flux.evaluate(solLeft, solRight), eNormal);
            }// for GP

            maxUL[iEdge] = uMaxOneEdge * e->length();
        }// for real edges 

       // collect result in each cell and compute tau
        for (size_t i = 0; i < M.nRealCells; ++i)
        {
            const shared_ptr<Cell> cell = M.cells[i];

            double uSum = 0.0;

            for (const shared_ptr<Edge> edge : cell->edges)
                uSum += maxUL[edge->number];

            factCo = tauOld * uSum / cell->getArea();
            relTau = CoNum / (factCo + 1e-6); //1e-6 is technical small number
            relTau = min(relTau, maxTauGrowth);
            tauNew = min(relTau * tauOld, maxTau);
            newTauLocal.push_back(tauNew);
        }

        tauNewLocal = *min_element(newTauLocal.begin(),newTauLocal.end());

        // collect with MPI to get min tau in full flow domain
        double tauNew = 0.0;
        MPI_Reduce(
            &tauNewLocal,
            &tauNew,
            1,
            MPI_DOUBLE,
            MPI_MIN,
            0,
            MPI_COMM_WORLD
        );

        MPI_Bcast(
            &tauNew,
            1,
            MPI_DOUBLE,
            0,
            MPI_COMM_WORLD
        );

        tauOld = tauNew;
        */
    
    }// if isDynamic

    t = tOld + tau;
    tOld = t;

}

bool TimeControl::isOutput()
{
    if (fabs(outputTime - t) < 1e-12)
    {
        timeListing << t << endl;
        outputTime += outputInterval;
        return true;
    }
    else
    {
        return false;
    }
}

