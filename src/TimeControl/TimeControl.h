#ifndef TIMECONTROL_H
#define TIMECONTROL_H

#include <fstream>

#include "Params.h"
#include "Mesh.h"
#include "Physics.h"
#include "Solution.h"

///
/// Store of time value, time step and controller of time step modification
///

class TimeControl
{

private:

    //------ Time itself
    
    /// Previous time point
    double tOld;

    /// Current time
    double t;

    /// End time
    double tEnd;

    /// Output time
    double outputTime;

    /// Output interval
    double outputInterval;

    //----- Time steps

    /// Old time step -- with which we're working during the time integrating
    double tau;

    /// New time step -- for the next iteration
    double tauNew;

    /// Saved time step which 
    double tauSaved;

    /// True if it is small time step for output
    bool stepForOutput;
    
    //----- Time step renewal

    /// Courant number
    double CoNum;

    /// Max time step
    double maxTau;

    /// Max time growth factor
    double maxTauGrowth;

    /// Switch static/dynamic time step
    bool isDynamic;

    /// Reference to mesh
    const Mesh& M;

    // Reference to physics
    const Physics& physics;

    //Reference to solution
    const Solution& sln;

    /// Stream to time listing file
    //std::ofstream timeListing;

public:

    /// Constructor
    TimeControl(
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
        const std::string listingPath = "alphaCoeffs/times"
    ); 

    /// Destructor
    ~TimeControl();

    /// Get time
    double getTime() const { return t; }
    
    /// Set time
    void updateTime(double t_new) { t = t_new; }

    /// Get time step
    double getTau() const { return tau; }

    /// Get new time step
    double getNewTau() const {return tauNew;}

    /// Update time step
    void updateTimeStep();

    /// Update current time in Runge --- Kutta step
    void updateTimeValueRK(double fracTau);

    /// Check if time is not finished
    bool running() { return t < tEnd; }

    /// Check time point for output
    bool isOutput();

};

#endif // TIMECONTROL_H
