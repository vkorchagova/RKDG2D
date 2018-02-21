#ifndef TIME_H
#define TIME_H


class Time
{
    //- Current time
    double runTime_;


public:

    //- Default constructor
    Time(){}

    //- Set time
    void updateTime (double t) {runTime_ = t;}

    //- Get time
    double runTime() const {return runTime_;}
};

#endif // TIME_H
