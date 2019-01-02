#ifndef TIMECONTROL_H
#define TIMECONTROL_H

#include "Params.h"
#include "Mesh.h"


class TimeControl
{
private:

	/// Time itself
	//- Current time
	double t;

	/// Time steps
    //- Old time step -- with which we're working during the time integrating
    double tau;
    //- New time step -- for the next iteration
    double tauNew;
	
	/// Time step renewal
    //- Courant number
    double CoNum;
	//- Max solution speed estimation
	double MaxSpeed;
    //- Max time step
    double maxTau;
    //- Max time growth factor
    double maxTauGrowth;
    //- Switch static/dynamic time step
    bool isDynamic;

	/// References to outer classes
    //- Reference to mesh
    const Mesh& M;

public:

    //- Constructor
    TimeControl(const Mesh& msh) : M(msh)
    {
		t=0.0;
        tau=1e-4;

        tauNew = tau;
		//write the construction with initialization from Params.h or some other data source;
	}  


	//- Get time
	double getTime() const { return t; }
	//- Set time
	void updateTime(double t_new) { t = t_new; }

	//- Get time step
	double getTau() const { return tau; }
    //- Get new time step
    double getNewTau() const {return tauNew;}
    //- Update time step
    void updateTimeStep(double MSpeed);

};

#endif // TIMECONTROL_H
