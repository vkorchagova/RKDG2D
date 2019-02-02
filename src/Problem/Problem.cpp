#include "Problem.h"
#include <iostream>
#include <omp.h>

using namespace std;

// ------------------ Constructors & Destructors ----------------

Problem::Problem(CaseInit task, const Mesh& m, const TimeControl& t) : M(m), T(t)
{
    setInitialConditions(task);
    setBoundaryConditions(task);
} // end constructor 

Problem::~Problem()
{

}

// ------------------ Private class methods --------------------

// ------------------ Public class methods --------------------

void Problem::setInitialConditions(CaseInit task)
{
	// Function for initial condition
	//double rho0 = 0.5;

	function<double(const Point& r)> initRho;
	function<double(const Point& r)> initP;
	function<double(const Point& r)> initU;
	function<double(const Point& r)> initV;


	switch (task)
	{

	case SodX:
	{
        cpcv = 1.4;

		initRho = [](const Point& r) { return (r.x() < 0) ? 1.0 : 0.125; };
		initP   = [](const Point& r) { return (r.x() < 0) ? 1.0 : 0.1;  };
		initU   = [](const Point& r) { return (r.x() < 0) ? 0.0 : 0.0; };
		initV   = [](const Point& r) { return 0.0; };

        break;
	}// SodX

	case SodY:
	{
        cpcv = 1.4;

		initRho = [](const Point& r) { return (r.y() < 0) ? 1.0 : 0.125; };
		initP = [](const Point& r) { return (r.y() < 0) ? 1.0 : 0.1;  };
		initV = [](const Point& r) { return (r.y() < 0) ? 0.75 : 0.0; };
		initU = [](const Point& r) { return 0.0; };

        break;
	}// SodY

	case SodDiag:
	{
        cpcv = 1.4;

		initRho = [](const Point& r) { return ((r.x() + r.y()) < 1.0) ? 1.0 : 0.125; };
		initP = [](const Point& r) { return ((r.x() + r.y()) < 1.0) ? 1.0 : 0.1;  };
		initV = [](const Point& r) { return 0.0; };
		initU = [](const Point& r) { return 0.0; };

        break;
	}// SodDiag

	case SodCircle:
	{
        cpcv = 1.4;

		initRho = [](const Point& r) { return (r.length() <= 0.1) ? 1.0 : 0.125; };
		initP = [](const Point& r) { return (r.length() <= 0.1) ? 1.0 : 0.1;  };
		initV = [](const Point& r) { return 0.0; };
		initU = [](const Point& r) { return 0.0; };

        break;
	}// SodCircle

    case ForwardStep:
    {
        cpcv = 1.4;

        initRho = [](const Point& r) { return 1.0; };
        initP   = [&](const Point& r) { return 1.0 / cpcv;  };
        initU   = [](const Point& r) { return 3.0; };
        initV   = [](const Point& r) { return 0.0; };

        break;
    } // forward step

    case Sedov:
    {
        cpcv = 1.4;

        initRho = [](const Point& r) { return 1.0; };
        initP   = [&](const Point& r) { return 1.0 / cpcv;  };
        initU   = [](const Point& r) { return 3.0; };
        initV   = [](const Point& r) { return 0.0; };

        break;
    } // forward step

	}// switch

    init = [=](const Point& r)
        {
            return numvector<double, 5> \
            { \
                initRho(r), \
                initRho(r)*initU(r), \
                initRho(r)*initV(r), \
                0.0, \
                initP(r) / (cpcv - 1.0) + \
                0.5 * initRho(r) * (sqr(initU(r)) + sqr(initV(r))) \
            };
        };

};// SetInitCond

	/*
    
    else if (caseName == "Woodward")
    {
        cpcv = 1.4;

        initRho = [](const Point& r) { return 1.0; };
        initP   = [](const Point& r) { return (r.x() <= -0.4) ? 1000.0 : ((r.x() >= 0.4) ? 100.0 : 0.01);  };
        initV   = [](const Point& r) { return 0.0; };
        initU   = [](const Point& r) { return 0.0; };
    }
    else if (caseName == "BlastCircle")
    {
	cpcv = 1.4;
    
	initRho = [](const Point& r) { return 1.0; };
	initP = [](const Point& r) { return (r.length() < 0.1) ? 10.0 : 0.1; };
	initU = [](const Point& r) { return 0.0; };
	initV = [](const Point& r) { return 0.0; };
    }
    else if (caseName == "Noh")
    {
        cpcv = 5.0/3.0;

        initRho = [](const Point& r) { return 1.0; };
        initP   = [](const Point& r) { return 1e-6;  };
        initU   = [](const Point& r) { return (r.x() < 0) ? 1.0 : -1.0; };
        initV   = [](const Point& r) { return 0.0; };
    }
    else if (caseName == "123")
    {
        cpcv = 1.4;

        initRho = [](const Point& r) { return 1.0; };
        initP   = [](const Point& r) { return 0.4;  };
        initU   = [](const Point& r) { return (r.x() < 0) ? -2.0 : 2.0; };
        initV   = [](const Point& r) { return 0.0; };
    }
    else if (caseName == "forwardStep")
    {
        cpcv = 1.4;

        initRho = [](const Point& r) { return 1.0; };
        initP   = [&](const Point& r) { return 1.0 / cpcv;  };
        initU   = [](const Point& r) { return 3.0; };
        initV   = [](const Point& r) { return 0.0; };
    }
    else if (caseName == "doubleMach")
    {
        cpcv = 1.4;

        initRho = [](const Point& r) { return (r.x() < 0.15) ? 8.0 : 1.4; };
        initP   = [](const Point& r) { return (r.x() < 0.15) ? 116.518 : 1.0;  };
        initU   = [](const Point& r) { return (r.x() < 0.15) ? 8.25 : 0.0; };
        initV   = [](const Point& r) { return 0.0; };
    }
    else if (caseName == "acousticPulse")
    {
        cpcv = 1.4;

        initRho = [](const Point& r) { return 1.0 + 1e-6*exp( - 40.0*sqr(r.x() )- 40.0*sqr(r.y() )); };
        initP   = [=](const Point& r) { return pow(initRho(r), cpcv);  };
        initU   = [](const Point& r) { return 0.0; };
        initV   = [](const Point& r) { return 0.0; };
    }
    else if (caseName == "monopole" || caseName == "dipole")
    {
        cpcv = 1.4;

        initRho = [](const Point& r) { return 1.0; };
        initP   = [&](const Point& r) { return 1.0 / cpcv;  };
        initU   = [](const Point& r) { return 0.0; };
        initV   = [](const Point& r) { return 0.0; };
    }
    else
    {
        cpcv = 1.4;

        cout << "Problem " << caseName << " not found\n";
        exit(0);
    }

    init = [=](const Point& r)
    {
        return numvector<double, 5> \
        { \
            initRho(r), \
            initRho(r)*initU(r), \
            initRho(r)*initV(r), \
            0.0, \
            initP(r) / (cpcv - 1.0) + \
                0.5 * initRho(r) * (sqr(initU(r)) + sqr(initV(r))) \
        };
    };

}
*/
void Problem::setBoundaryConditions(CaseInit task)
{
    // shared_ptr<BoundarySine> bSine = make_shared<BoundarySine>(1e-3,0.5,time,*this);

    //vector<shared_ptr<Boundary>> bc = {};
    switch (task)
    {

    case SodX:
    case SodY:
    case SodDiag:
    case SodCircle:
    case Sedov:
    {
        for (const Patch& p : M.patches)
            bc.emplace_back(make_shared<BoundarySlip>(p));

        break;
    }

    case ForwardStep:
    {
        bc.emplace_back(
            make_shared<BoundaryConstant>(
                M.patches[0], 
                numvector<double, dimPh>({1.0, 3.0, 0.0, 0.0, 6.286})
            )
        );
        bc.emplace_back(make_shared<BoundaryOpen>(M.patches[1]));
        bc.emplace_back(make_shared<BoundarySlip>(M.patches[2]));
        bc.emplace_back(make_shared<BoundarySlip>(M.patches[3]));

        break;
    }
    }
   
    
}


