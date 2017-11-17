#pragma once

#include "Problem.h"

class Flux
{
public:
	Flux();
	~Flux();

	virtual vector<double> evaluateHor() = 0;
	virtual vector<double> evaluateVer() = 0;
};

