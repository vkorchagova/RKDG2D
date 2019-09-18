#include "Problem.h"
#include <iostream>
#include <omp.h>

using namespace std;

// ------------------ Constructors & Destructors ----------------

Problem::Problem(CaseInit task, const Mesh& m, const TimeControl& t, Physics& phs) : M(m), T(t), phs(phs)
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

    source = [=](const numvector<double, dimPh> sol, const Point& r)
        {
            return numvector<double, 5> \
            { \
                0.0, \
                0.0, \
                0.0, \
                0.0, \
                0.0 \
            }; 
        }; 


	switch (task)
	{

	case SodX:
	{
        phs.cpcv = 1.4;
        phs.covolume = 0.0;

		initRho = [](const Point& r) { return (r.x() < 0) ? 1.0 : 0.125; };
		initP   = [](const Point& r) { return (r.x() < 0) ? 1.0 : 0.1;  };
		initU   = [](const Point& r) { return (r.x() < 0) ? 0.0 : 0.0; };
		initV   = [](const Point& r) { return 0.0; };

        break;
	}// SodX

    case SodXCovol:
    {
        phs.cpcv = 1.3;
        phs.covolume = 0.001;

        initRho = [](const Point& r) { return (r.x() < -0.1) ? 100.0 : 1.0; };
        initP   = [](const Point& r) { return (r.x() < -0.1) ? 100.0 * 1e6 : 0.1 * 1e6;  };
        initU   = [](const Point& r) { return (r.x() < -0.1) ? 0.0 : 0.0; };
        initV   = [](const Point& r) { return 0.0; };

        break;
    }// SodX

	case SodY:
	{
        phs.cpcv = 1.4;
        phs.covolume = 0.0;

		initRho = [](const Point& r) { return (r.y() < 0) ? 1.0 : 0.125; };
		initP = [](const Point& r) { return (r.y() < 0) ? 1.0 : 0.1;  };
		initV = [](const Point& r) { return (r.y() < 0) ? 0.75 : 0.0; };
		initU = [](const Point& r) { return 0.0; };

        break;
	}// SodY

	case SodDiag:
	{
        phs.cpcv = 1.4;
        phs.covolume = 0.0;

		initRho = [](const Point& r) { return ((r.x() + r.y()) < 1.0) ? 1.0 : 0.125; };
		initP = [](const Point& r) { return ((r.x() + r.y()) < 1.0) ? 1.0 : 0.1;  };
		initV = [](const Point& r) { return 0.0; };
		initU = [](const Point& r) { return 0.0; };

        break;
	}// SodDiag

	case SodCircle:
	{
        phs.cpcv = 1.4;
        phs.covolume = 0.0;

		initRho = [](const Point& r) { return (r.length() <= 0.2) ? 1.0 : 0.125; };
		initP = [](const Point& r) { return (r.length() <= 0.2) ? 1.0 : 0.1;  };
		initV = [](const Point& r) { return 0.0; };
		initU = [](const Point& r) { return 0.0; };

        break;
	}// SodCircle

    case ForwardStep:
    {
        phs.cpcv = 1.4;
        phs.covolume = 0.0;

        initRho = [](const Point& r) { return 1.0; };
        initP   = [&](const Point& r) { return 1.0 / phs.cpcv;  };
        initU   = [](const Point& r) { return 3.0; };
        initV   = [](const Point& r) { return 0.0; };

        break;
    } // forward step

    case Sedov:
    {
        phs.cpcv = 1.4;
        phs.covolume = 0.0;

        initRho = [](const Point& r) { return 1.0; };
        initP   = [&](const Point& r) { return 1.0 / phs.cpcv;  };
        initU   = [](const Point& r) { return 3.0; };
        initV   = [](const Point& r) { return 0.0; };

        break;
    } // Sedov

    case DoubleMach:
    {
        phs.cpcv = 1.4;
        phs.covolume = 0.0;

        initRho = [](const Point& r) { return (r.x() < 0.15) ? 8.0 : 1.4; };
        initP   = [](const Point& r) { return (r.x() < 0.15) ? 116.518 : 1.0;  };
        initU   = [](const Point& r) { return (r.x() < 0.15) ? 8.25 : 0.0; };
        initV   = [](const Point& r) { return 0.0; };

        break;
    } // DMach

    case Ladenburg:
    {
        phs.cpcv = 1.4;
        phs.covolume = 0.0;

        initRho = [](const Point& r) { return 1.184308158; };
        initP   = [](const Point& r) { return 101325.0;  };
        initU   = [](const Point& r) { return 0.0; };
        initV   = [](const Point& r) { return 0.0; };

        break;
    } // Ladenburg

    case ShuOsher:
    {
        phs.cpcv = 1.4;
        phs.covolume = 0.0;

        initRho = [](const Point& r) { return (r.x() < -0.5 + 0.125) ? 3.857143 : 1.0 + 0.2 * sin(8.0 * 2.0 * 3.14159265*r.x() ); };
        initP   = [](const Point& r) { return (r.x() < -0.5 + 0.125) ? 10.3333 : 1.0;  };
        initU   = [](const Point& r) { return (r.x() < -0.5 + 0.125) ? 2.629369 : 0.0; };
        initV   = [](const Point& r) { return 0.0; };

        break;
    }

	}// switch

    init = [=](const Point& r)
        {
            return numvector<double, 5> \
            { \
                initRho(r), \
                initRho(r)*initU(r), \
                initRho(r)*initV(r), \
                0.0, \
                phs.e(initRho(r), initU(r), initV(r), 0.0, initP(r)) \
            }; 
        }; //for Co-Volume 0.001, for ideal 0.0


};// SetInitCond

	/*
    
    else if (caseName == "Woodward")
    {
        phs.cpcv = 1.4;

        initRho = [](const Point& r) { return 1.0; };
        initP   = [](const Point& r) { return (r.x() <= -0.4) ? 1000.0 : ((r.x() >= 0.4) ? 100.0 : 0.01);  };
        initV   = [](const Point& r) { return 0.0; };
        initU   = [](const Point& r) { return 0.0; };
    }
    else if (caseName == "BlastCircle")
    {
	phs.cpcv = 1.4;
    
	initRho = [](const Point& r) { return 1.0; };
	initP = [](const Point& r) { return (r.length() < 0.1) ? 10.0 : 0.1; };
	initU = [](const Point& r) { return 0.0; };
	initV = [](const Point& r) { return 0.0; };
    }
    else if (caseName == "Noh")
    {
        phs.cpcv = 5.0/3.0;

        initRho = [](const Point& r) { return 1.0; };
        initP   = [](const Point& r) { return 1e-6;  };
        initU   = [](const Point& r) { return (r.x() < 0) ? 1.0 : -1.0; };
        initV   = [](const Point& r) { return 0.0; };
    }
    else if (caseName == "123")
    {
        phs.cpcv = 1.4;

        initRho = [](const Point& r) { return 1.0; };
        initP   = [](const Point& r) { return 0.4;  };
        initU   = [](const Point& r) { return (r.x() < 0) ? -2.0 : 2.0; };
        initV   = [](const Point& r) { return 0.0; };
    }
    else if (caseName == "forwardStep")
    {
        phs.cpcv = 1.4;

        initRho = [](const Point& r) { return 1.0; };
        initP   = [&](const Point& r) { return 1.0 / phs.cpcv;  };
        initU   = [](const Point& r) { return 3.0; };
        initV   = [](const Point& r) { return 0.0; };
    }
    else if (caseName == "doubleMach")
    {
        phs.cpcv = 1.4;

        initRho = [](const Point& r) { return (r.x() < 0.15) ? 8.0 : 1.4; };
        initP   = [](const Point& r) { return (r.x() < 0.15) ? 116.518 : 1.0;  };
        initU   = [](const Point& r) { return (r.x() < 0.15) ? 8.25 : 0.0; };
        initV   = [](const Point& r) { return 0.0; };
    }
    else if (caseName == "acousticPulse")
    {
        phs.cpcv = 1.4;

        initRho = [](const Point& r) { return 1.0 + 1e-6*exp( - 40.0*sqr(r.x() )- 40.0*sqr(r.y() )); };
        initP   = [=](const Point& r) { return pow(initRho(r), phs.cpcv);  };
        initU   = [](const Point& r) { return 0.0; };
        initV   = [](const Point& r) { return 0.0; };
    }
    else if (caseName == "monopole" || caseName == "dipole")
    {
        phs.cpcv = 1.4;

        initRho = [](const Point& r) { return 1.0; };
        initP   = [&](const Point& r) { return 1.0 / phs.cpcv;  };
        initU   = [](const Point& r) { return 0.0; };
        initV   = [](const Point& r) { return 0.0; };
    }
    else
    {
        phs.cpcv = 1.4;

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
            initP(r) / (phs.cpcv - 1.0) + \
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
    case SodXCovol:
    case SodY:
    case SodDiag:
    case SodCircle:
    case Sedov:
    {
        for (const Patch& p : M.patches)
            bc.emplace_back(make_shared<BoundarySlip>(p));

        break;
    }
    
    case ShuOsher:
    {
        for (const Patch& p : M.patches)
            bc.emplace_back(make_shared<BoundaryOpen>(p));

        break;
    }

    case ForwardStep:
    {
        double inletRho = 1.0;
        double inletU = 3.0;
        double inletV = 0.0;
        double inletP = 1.0 / phs.cpcv;

        for (const Patch& patch : M.patches)
        {
            if (patch.name == "left" or 
                patch.name == "inlet")
            {
                bc.emplace_back(
                make_shared<BoundaryConstant>(
                    patch, 
                    numvector<double, dimPh>(
                    {
                        inletRho, 
                        -inletU * inletRho, 
                         inletV * inletRho, 
                        0.0, 
                        phs.e(inletRho, inletU, inletV, 0.0, inletP)
                    })
                ));
            }
            else if (patch.name == "right" or 
                     patch.name == "outlet")
            {
                bc.emplace_back(make_shared<BoundaryOpen>(patch));
            }
            else if (patch.name == "up" or 
                     patch.name == "top" or 
                     patch.name == "down" or
                     patch.name == "walls" or 
                     patch.name == "wall" or  
                     patch.name == "step")
            {
                bc.emplace_back(make_shared<BoundarySlip>(patch));
            }
            else
            {
                cout << "Boundary condition for patch " << patch.name << " is not found.\n"
                     << "Check settings in src/Problem.cpp/setBoundaryConditions" << endl;
                exit(1);
            }

            break;
        }
    }

    case DoubleMach:
    {
        double inletRho = 8.0;
        double inletU = 8.25;
        double inletV = 0.0;
        double inletP = 116.518;

        for (const Patch& patch : M.patches)
        {
            if (patch.name == "left" or 
                patch.name == "inlet")
            {
                bc.emplace_back(
                make_shared<BoundaryConstant>(
                    patch, 
                    numvector<double, dimPh>(
                    {
                        inletRho, 
                        -inletU * inletRho, 
                         inletV * inletRho, 
                        0.0, 
                        phs.e(inletRho, inletU, inletV, 0.0, inletP)
                    })
                ));
            }
            else if (patch.name == "right" or 
                     patch.name == "outlet")
            {
                bc.emplace_back(make_shared<BoundaryOpen>(patch));
            }
            else if (patch.name == "up" or 
                     patch.name == "top" or 
                     patch.name == "walls" or 
                     patch.name == "wall" or 
                     patch.name == "down")
            {
                bc.emplace_back(make_shared<BoundarySlip>(patch));
            }
            else
            {
                cout << "Boundary condition for patch " << patch.name << " is not found.\n"
                     << "Check settings in src/Problem.cpp/setBoundaryConditions" << endl;
                exit(1);
            }
        }
        break;
    }

    case Ladenburg:
    {
        double inletRho = 3.8303970021; //3.830183898
        double inletU = 315.6;
        double inletV = 0.0;
        double inletP = 271724;

        for (const Patch& patch : M.patches)
        {
            if (patch.name == "left" or 
                patch.name == "inlet")
            {
                bc.emplace_back(
                make_shared<BoundaryConstant>(
                    patch, 
                    numvector<double, dimPh>(
                    {
                        inletRho, 
                        -inletU * inletRho, 
                         inletV * inletRho, 
                        0.0, 
                        phs.e(inletRho, inletU, inletV, 0.0, inletP)
                    })
                ));
            }
            else if (patch.name == "right" or 
                     patch.name == "outlet")
            {
                bc.emplace_back(make_shared<BoundaryOpen>(patch));
            }
            else if (patch.name == "up" or 
                     patch.name == "top" or 
                     patch.name == "down" or
                     patch.name == "walls" or 
                     patch.name == "wall")
            {
                bc.emplace_back(make_shared<BoundarySlip>(patch));
            }
            else
            {
                cout << "Boundary condition for patch " << patch.name << " is not found.\n"
                     << "Check settings in src/Problem.cpp/setBoundaryConditions" << endl;
                exit(1);
            }
        }
        
        break;
    }
    
    }
}


