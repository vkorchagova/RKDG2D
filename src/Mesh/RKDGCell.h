#ifndef RKDGCELL_H
#define RKDGCELL_H


#include "Cell.h"
#include "Problem.h"

class RKDGCell : public Cell
{

public:

    //- Pointer to problem
    Problem* problem;


public:

    //- Default constructor
    RKDGCell();

    //- Destructor
    ~RKDGCell();

    //- Calculate local RHS
    double getLocalRHS();


}; // RKDGCell

#endif // RKDGCELL_H
