#ifndef TIMECLASS_H
#define TIMECLASS_H

/// Time object
///

class Time
{
    /// Current time
    double runTime_;
    
    /// Current time step
    double tau;


public:

    /// Default constructor
    Time(){}

    /// Set time
    void updateTime (double t) {runTime_ = t;}

    /// Get time
    double runTime() const {return runTime_;}
};

#endif // TIME_H
