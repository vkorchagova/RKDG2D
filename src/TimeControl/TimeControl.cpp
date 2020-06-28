#include "TimeControl.h"
#include <algorithm>
#include <iomanip>

using namespace std;



TimeControl::TimeControl(
        const Mesh& msh,
        const Physics& phs,
        const Solution& sln, 
        double tStart, 
        double tEnd, 
        double initTau, 
        double outputInterval,
        bool   isDynamic,
        double CoNum,
        double maxTau,
        double maxTauGrowth,
        const string listingPath
    ) : 
    M (msh), 
    physics (phs),
    sln (sln),
    t (tStart), 
    tau (initTau), 
    tEnd (tEnd), 
    isDynamic (isDynamic),
    CoNum (CoNum),
    maxTau (maxTau),
    maxTauGrowth (maxTauGrowth),
    outputInterval (outputInterval)
{
    tauNew = tau;
    tauSaved = initTau;//0.0;
    tOld = tStart;

    outputTime = tStart + outputInterval;
    //timeListing.open(listingPath.c_str());

    //if (!timeListing.is_open())
    //{
    //    cout << "File " << listingPath << " could not be opened\n";
    //    exit(0);
    //}

    //timeListing << setprecision(6);
    //timeListing << 0.000000 << endl;
    //write the construction with initialization from Params.h or some other data source;
}

TimeControl::~TimeControl()
{
    //timeListing.close();
}

void TimeControl::updateTimeValueRK(double fracTau)
{
    t = tOld + fracTau;
}

void TimeControl::updateTimeStep() 
{
    if (isDynamic) // template for the code of Vik
    {
        double factCo = 0.0;
        double relTau = 0.0;
        double tauNewLocal = 0.0;
        double tauNewCell = 0.0;

        vector<double> newTauLocal;
        newTauLocal.reserve(M.nRealCells);

        vector<double> maxU;
        maxU.reserve(M.nRealCells);

        // collect result in each cell and compute tau
        for (size_t i = 0; i < M.nRealCells; ++i)
        {
            const shared_ptr<Cell> cell = M.cells[i];

            double uMax = 0.0;

            for (const shared_ptr<Point> p : cell->nodes)
            {
                numvector<double, dimPh> solNode = sln.reconstruct(cell->number,*p);
                uMax = max(uMax, physics.magU(solNode));
            }

            double perimeter = 0.0;

            for (const shared_ptr<Edge> edge : cell->edges)
                perimeter += edge->getLength();

            factCo = tau * perimeter * uMax / cell->getArea();
            relTau = CoNum / (factCo + 1e-6); //1e-6 is technical small number
            relTau = max(relTau, 1.0 - maxTauGrowth);
            relTau = min(relTau, 1.0 + maxTauGrowth);
            tauNewCell = min(relTau * tau, maxTau);
            newTauLocal.push_back(tauNewCell);
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

        if (myRank == 0)
        {
            tauNew = min(tauNew, outputTime - t);

            if (tauNew < (1.0 - maxTauGrowth) * tau && !stepForOutput)
            {
                tauSaved = tau;
                stepForOutput = true;
            }
        }

        MPI_Bcast(
            &tauNew,
            1,
            MPI_DOUBLE,
            0,
            MPI_COMM_WORLD
        );

        tau = tauNew;
    
    }// if isDynamic

    tOld = t;
}

bool TimeControl::isOutput()
{
    if (fabs(outputTime - t) < 1e-7)
    {
        //timeListing << t << endl;
        outputTime += outputInterval;

        if (stepForOutput)
        {
            stepForOutput = false;
            tau = tauSaved;
        }
        return true;
    }
    else
    {
        return false;
    }
}

