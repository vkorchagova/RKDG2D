#ifndef TIMECONTROL_H
#define TIMECONTROL_H

#include "Mesh2D.h"

/// Dynamic time controller
///
/// Change time step value according to CFL number


class TimeControl
{
private:

    /// Old time step
    double tauOld;

    /// New time step
    double tauNew;

    /// Constant reference to mesh
    const Mesh2D& mesh;

    /// Courant number
    double CoNum;

    /// Max time step
    double maxTau;

    /// Max time growth factor
    double maxTauGrowth;

    /// Switch static/dynamic time step
    bool modifyTime;

    /// Mass flux in all studies of solution
    std::vector<double> massFlux;

public:

    /// Constructor
    TimeControl(const Mesh2D& msh, double Co, double maxDeltaT, double maxGrowth, double initTau, bool modified = false) : \
                mesh(msh), CoNum(Co), maxTau(maxDeltaT), maxTauGrowth(maxGrowth), tauOld(initTau), modifyTime(modified) \
                {tauNew = tauOld; massFlux.resize(mesh.nCells); std::fill(massFlux.begin(), massFlux.end(), 0.0);}

    /// Get new time step
    double getNewTau() const {return tauNew;}

    /// Update time step
    void updateTimeStep();

    /// Get fabs of mass flux through cell edges in one study
    void getMassFlux();

};

#endif // TIMECONTROL_H
