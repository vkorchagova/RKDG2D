#include "Problem.h"
#include <iostream>
#include <omp.h>

using namespace std;

// ------------------ Constructors & Destructors ----------------

Problem::Problem(CaseInit task, const Mesh& m, const TimeControl& t) : M(m), T(t)
{
    setInitialConditions(task);
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

void Problem::setBoundaryConditions(string caseName, const std::vector<Patch>& patches)
{
    // shared_ptr<BoundarySine> bSine = make_shared<BoundarySine>(1e-3,0.5,time,*this);

    vector<shared_ptr<Boundary>> bc = {};

    if (caseName == "SodX" || \
        caseName == "SodY" || \
        caseName == "SodDiag" || \
        caseName == "Woodward" || \
        caseName == "BlastCircle" || \
        caseName == "Noh" || \
        caseName == "123" || \
        caseName == "acousticPulse")
    {
        shared_ptr<BoundaryOpen> bOpen = make_shared<BoundaryOpen>();
        shared_ptr<BoundarySlip> bSlip = make_shared<BoundarySlip>();

        //bc = {bOpen, bOpen};

        bc = {bSlip, bSlip, bSlip, bSlip};

        //bc = {bOpen, bOpen, bOpen, bOpen};
    }
    else if (caseName == "forwardStep")
    {
        shared_ptr<BoundaryOpen> bOpen = make_shared<BoundaryOpen>();
        shared_ptr<BoundarySlip> bSlip = make_shared<BoundarySlip>();
        shared_ptr<BoundaryConstant> bConst = \
                make_shared<BoundaryConstant>( \
                numvector<double, 5>({1.0, -3.0, 0.0, 0.0, 6.286}));

        bc = {bConst, bOpen, bSlip, bSlip};
    }    
    else if (caseName == "doubleMach")
    {
        shared_ptr<BoundaryOpen> bOpen = make_shared<BoundaryOpen>();
        shared_ptr<BoundarySlip> bSlip = make_shared<BoundarySlip>();
        shared_ptr<BoundaryConstant> bConst = 
                make_shared<BoundaryConstant>( \
                numvector<double, 5>({8.0, -8.25*8.0, 0.0, 0.0, 563.545}));

        bc = {bConst, bOpen, bSlip};
    }
    else if (caseName == "SodCircle")
    {
        shared_ptr<BoundarySlip> bSlip = make_shared<BoundarySlip>();

        bc = {bSlip, bSlip, bSlip};
    }
    else if (caseName == "monopole")
    {
        //shared_ptr<BoundaryOpen> bOpen = make_shared<BoundaryOpen>();
        //shared_ptr<BoundarySlip> bSlip = make_shared<BoundarySlip>();
        shared_ptr<BoundarySine> bSine = \
                make_shared<BoundarySine>(1e-6,5.0,time,*this,init(Point({0.0,0.0})));
        shared_ptr<BoundaryConstant> bConst = \
                make_shared<BoundaryConstant>(init(Point({0.0,0.0})));


        bc = {bSine, bConst};
    }
    else if (caseName == "dipole")
    {
        //shared_ptr<BoundaryOpen> bOpen = make_shared<BoundaryOpen>();
        //shared_ptr<BoundarySlip> bSlip = make_shared<BoundarySlip>();
        shared_ptr<BoundarySineDir> bSine = \
                make_shared<BoundarySineDir>(1e-6,5.0,time,*this,init(Point({0.0,0.0})));
        shared_ptr<BoundaryConstant> bConst = \
                make_shared<BoundaryConstant>(init(Point({0.0,0.0})));


        bc = {bSine, bConst};
    }
    else
    {
        cout << "Problem " << caseName << " not found\n";
        exit(0);
    }


    for (int i = 0; i < patches.size(); ++i)
        for (int j = 0; j < patches[i].edgeGroup.size(); ++j)
            patches[i].edgeGroup[j]->setBoundary(bc[i]);

    for (int i = 0; i < patches.size(); ++i)
        cout << "Patch #" << i << ": type = " << bc[i]->type << endl;
}
*/

