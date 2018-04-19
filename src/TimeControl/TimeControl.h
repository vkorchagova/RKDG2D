#ifndef TIMECONTROL_H
#define TIMECONTROL_H

#include "Mesh2D.h"


class TimeControl
{
private:

    //- Old time step
    double tauOld;

    //- New time step
    double tauNew;

    //- Reference to mesh
    const Mesh2D& mesh;

    //- Courant number
    double CoNum;

    //- Max time step
    double maxTau;

    //- Max time growth factor
    double maxTauGrowth;

    //- Switch static/dynamic time step
    bool modifyTime;

public:

    //- Constructor
    TimeControl(const Mesh2D& msh, double Co, double maxDeltaT, double maxGrowth, double initTau, bool modified = false) : \
                mesh(msh), CoNum(Co), maxTau(maxDeltaT), maxTauGrowth(maxGrowth), tauOld(initTau), modifyTime(modified) \
                {tauNew = tauOld;}

    //- Get new time step
    double getNewTau() const {return tauNew;}

    //- Update time step
    void updateTimeStep();

};

#endif // TIMECONTROL_H
