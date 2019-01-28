#include "TimeControl.h"
#include <algorithm>
#include <iomanip>

using namespace std;



TimeControl::TimeControl(
        const Mesh& msh, 
        const double tStart, 
        const double tEnd, 
        const double initTau, 
        const double outputInterval,
        const string listingPath
    ) : 
    M(msh), 
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
    /*
    if (isDynamic) // template for the code of Vik
    {
        double factCo = 0.0;
        double relTau = 0.0;

        vector<double> newTauLocal;
        newTauLocal.reserve(M.nRealCells);

		MaxSpeed = MSpeed;

//#pragma omp parallel for \
shared (M, newTauLocal) \
private (factCo, relTau) \
default (none)
        //for (const shared_ptr<Cell> cell : mesh.cells)

		??? CYCLE THROUGH THE CELLS IS NEEDED???YESSSS 
        factCo = tau * MaxSpeed / cell->getArea();
        relTau = CoNum / (factCo + 1e-6); //1e-6 is technical small number
        relTau = min(relTau, maxTauGrowth);
        tauNew = min(relTau * tau, maxTau);
        ???newTauLocal.push_back(tauNew);

        tauNew = *min_element(newTauLocal.begin(),newTauLocal.end());

        tau = tauNew;
    }// if isDynamic

	if (isDynamic)
	{
		///collection
		// Firstly we need the mesh "size"
		double H = 0.0; // !!! if this way, need it from the mesh constructor

		// Second -- the medium speed
		MaxSpeed=MSpeed;

		tau = CoNum*H / MaxSpeed;

	}// if isDynamic
*/
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

