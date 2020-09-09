#include "Problem.h"
#include <iostream>
#include <omp.h>

#include <cmath>
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
        phs.mu = 1e-3;
        phs.Pr = 1.0;

        initRho = [](const Point& r) { return (r.x() < 0) ? 1.0: 0.125; };
        initP = [](const Point& r) { return (r.x() < 0) ? 1.0 : 0.1;  };
        initU = [](const Point& r) { return (r.x() < 0) ? 0.0 : 0.0; };
		initV = [](const Point& r) { return 0.0; };

		break;
	}// SodX

	case SodXCovol:
	{
		phs.cpcv = 1.3;
		phs.covolume = 0.001;

		initRho = [](const Point& r) { return (r.x() < -0.1) ? 100.0 : 1.0; };
		initP = [](const Point& r) { return (r.x() < -0.1) ? 100.0 * 1e6 : 0.1 * 1e6;  };
		initU = [](const Point& r) { return (r.x() < -0.1) ? 0.0 : 0.0; };
		initV = [](const Point& r) { return 0.0; };

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
        phs.mu = 0.0;

        initRho = [](const Point& r) { return (r.length() <= 0.4) ? 1.0 : 0.125; };
        initP = [](const Point& r) { return (r.length() <= 0.4) ? 1.0 : 0.1;  };
		initV = [](const Point& r) { return 0.0; };
		initU = [](const Point& r) { return 0.0; };

		break;
	}// SodCircle

	case ForwardStep:
	{
		phs.cpcv = 1.4;
		phs.covolume = 0.0;

		initRho = [](const Point& r) { return 1.0; };
		initP = [&](const Point& r) { return 1.0 / phs.cpcv;  };
		initU = [](const Point& r) { return 3.0; };
		initV = [](const Point& r) { return 0.0; };

		break;
	} // forward step

	case Sedov:
	{
		phs.cpcv = 1.4;
		phs.covolume = 0.0;

		initRho = [](const Point& r) { return 1.0; };
		initP = [&](const Point& r) { return 1.0 / phs.cpcv;  };
		initU = [](const Point& r) { return 3.0; };
		initV = [](const Point& r) { return 0.0; };

		break;
	} // Sedov

	case DoubleMach:
	{
		phs.cpcv = 1.4;
		phs.covolume = 0.0;

		initRho = [](const Point& r) { return (r.x() < 0.15) ? 8.0 : 1.4; };
		initP = [](const Point& r) { return (r.x() < 0.15) ? 116.518 : 1.0;  };
		initU = [](const Point& r) { return (r.x() < 0.15) ? 8.25 : 0.0; };
		initV = [](const Point& r) { return 0.0; };

		break;
	} // DMach

	case Ladenburg:
	{
		phs.cpcv = 1.4;
		phs.covolume = 0.0;

		initRho = [](const Point& r) { return 1.184308158; };
		initP = [](const Point& r) { return 101325.0;  };
		initU = [](const Point& r) { return 0.0; };
		initV = [](const Point& r) { return 0.0; };

		break;
	} // Ladenburg

	case ShuOsher:
	{
		phs.cpcv = 1.4;
		phs.covolume = 0.0;

		initRho = [](const Point& r) { return (r.x() < -0.5 + 0.125) ? 3.857143 : 1.0 + 0.2 * sin(8.0 * 2.0 * 3.14159265*r.x()); };
		initP = [](const Point& r) { return (r.x() < -0.5 + 0.125) ? 10.3333 : 1.0;  };
		initU = [](const Point& r) { return (r.x() < -0.5 + 0.125) ? 2.629369 : 0.0; };
		initV = [](const Point& r) { return 0.0; };

		break;
	}

    case Blasius:
    {
        phs.cpcv = 1.4;
        phs.covolume = 0.0;
        phs.Pr = 100;//0.135;
        phs.mu = 1.25e-1;


        initRho = [](const Point& r) { return 1.4; };
        initP = [&](const Point& r) { return 1.0;  };
        initU = [](const Point& r) { return 0.1; };
        initV = [](const Point& r) { return 0.0; };

        break;
    } // Blasius problem

case Vortex:
	{
		phs.cpcv = 1.4;
		phs.covolume = 0.0;
		phs.Pr = 100;//0.135;
		phs.mu = 1.00e-2;
		phs.g = 1.0;




		initRho = [](const Point& r) { return 1.0; };
		initP = [=](const Point& r) {
			double mn, exp2, exp4, cft, ei2, ei4, res;

			const double myrho = 1.0;
			const double t0 = 0.2 * myrho / phs.mu;
			
			double qqq = r.length();

			mn = phs.g * phs.g * myrho / (8.0 * M_PI * M_PI * qqq * qqq);
			exp2 = std::exp(-qqq * qqq / (2.0 * phs.mu / myrho * t0));
			exp4 = std::exp(-qqq * qqq / (4.0 * phs.mu / myrho * t0));
			
			cft = phs.g * phs.g * myrho / (16.0 * M_PI * M_PI * phs.mu / myrho * t0);
			ei2 = std::expint(-qqq * qqq / (2.0 * phs.mu / myrho * t0));
			ei4 = std::expint(-qqq * qqq / (4.0 * phs.mu / myrho * t0));

			res = -mn - mn * exp2 + 2.0 * mn * exp4 - cft * ei2 + cft * ei4;
			
			return 1.0e+5 + res; };
						
		initU = [=](const Point& r) { 
			const double myrho = 1.0;
			const double t0 = 0.2 * myrho / phs.mu;
			double beta = std::atan2(r.y(), r.x());	

			double Vphi = 0.5 * phs.g / (M_PI * r.length()) * (1.0 - std::exp(-r.length()*r.length() / (4.0 * phs.mu / myrho * t0)));

			return (-Vphi * sin(beta)); };
		
		initV = [=](const Point& r) {
			const double myrho = 1.0;
			const double t0 = 0.2 * myrho / phs.mu;
			double beta = std::atan2(r.y(), r.x());

			double Vphi = 0.5 * phs.g / (M_PI * r.length()) * (1.0 - std::exp(-r.length()* r.length() / (4.0 * phs.mu / myrho * t0)));

			return (Vphi * cos(beta)); };

		break;
	} // vortex problem
	
case CylinderFlow:
    {
        phs.cpcv = 1.4;
        phs.covolume = 0.0;
        phs.Pr = 100;//0.135;
        phs.mu = 1.25e-1;


        initRho = [](const Point& r) { return 1.4; };
        initP = [&](const Point& r) { return 1.0;  };
        initU = [](const Point& r) { return 0.1; };
        initV = [](const Point& r) { return 0.0; };

        break;
    } // CylinderFlow problem


case AstroTest:
	{
		phs.cpcv = 5.0 / 3.0;
		phs.covolume = 0.0;

		double u = 1.02; // radial velocity
		double kTilde = 5.5e6;
		double M0 = 1e30;
		double G0 = 6.67e-11;
		double Ror = 696e6;

		double k = kTilde * pow(M0, phs.cpcv - 2.0) / G0 / pow(Ror, 3.0*phs.cpcv - 4.0);
		double frackgamma = (phs.cpcv - 1.0) / k / phs.cpcv;

		source = [=](const numvector<double, dimPh> sol, const Point& r)
		{
			return numvector<double, 5> \
			{ \
				0.0, \
				- r.x() / pow(r.length(), 3), \
				- r.y() / pow(r.length(), 3), \
				0.0, \
				0.0 \
			};
		};

		initRho = [=](const Point& r) { return     pow(frackgamma * (u*u - 1.0) * (r.length() - 1.0) / r.length(), 1.0 / (phs.cpcv - 1.0)); };
		initP = [=](const Point& r) { return k * pow(frackgamma * (u*u - 1.0) * (r.length() - 1.0) / r.length(), phs.cpcv / (phs.cpcv - 1.0));  };
		initU = [=](const Point& r) { return -u * r.y() / pow(r.length2(), 0.75); };
		initV = [=](const Point& r) { return   u * r.x() / pow(r.length2(), 0.75); };

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

	case AstroTest:
	{
		for (const Patch& p : M.patches)
			bc.emplace_back(make_shared<BoundarySlip>(p));
		/*
		double u = 1.02; // radial velocity
		double kTilde = 5.5e6;
		double M0 = 1e30;
		double G0 = 6.67e-11;
		double Ror = 696e6;

		double k = kTilde * pow(M0, phs.cpcv - 2.0) / G0 / pow(Ror, 3.0*phs.cpcv - 4.0);
		double frackgamma = (phs.cpcv - 1.0) / k / phs.cpcv;



		for (const Patch& patch : M.patches)
		{
			if (patch.name == "inner")
			{
				double inletRho = pow( frackgamma * (u*u - 1.0) * 0.1 / 1.1, 1.0 / (phs.cpcv - 1.0));
				double inletU = 0.0;
				double inletV = u / sqrt(1.1);
				double inletP = k * pow( frackgamma * (u*u - 1.0) * 0.1 / 1.1, phs.cpcv / (phs.cpcv - 1.0));

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
			else if (patch.name == "outer")
			{
				double inletRho = pow( frackgamma * (u*u - 1.0) * 2.1 / 3.1, 1.0 / (phs.cpcv - 1.0));
				double inletU = 0.0;
				double inletV = u / sqrt(3.1);
				double inletP = k * pow( frackgamma * (u*u - 1.0) * 2.1 / 3.1, phs.cpcv / (phs.cpcv - 1.0));

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
			else
			{
				cout << "Boundary condition for patch " << patch.name << " is not found.\n"
					 << "Check settings in src/Problem.cpp/setBoundaryConditions" << endl;
				exit(1);
			}
		}
		*/

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
			if (patch.name == "left" ||
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
			else if (patch.name == "right" ||
				patch.name == "outlet")
			{
				bc.emplace_back(make_shared<BoundaryOpen>(patch));
			}
			else if (patch.name == "up" ||
				patch.name == "top" ||
				patch.name == "down" ||
				patch.name == "walls" ||
				patch.name == "wall" ||
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
		}
		break;
	}

	case Vortex:
	{
		for (const Patch& patch : M.patches)
		{
			bc.emplace_back(make_shared<BoundarySlip>(patch));
		}
		break;
	}
	case DoubleMach:
	{
		double inletRho = 8.0;
		double inletU = 8.25;
		double inletV = 0.0;
		double inletP = 116.518;

		for (const Patch& patch : M.patches)
		{
			if (patch.name == "left" ||
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
			else if (patch.name == "right" ||
				patch.name == "outlet")
			{
				bc.emplace_back(make_shared<BoundaryOpen>(patch));
			}
			else if (patch.name == "up" ||
				patch.name == "top" ||
				patch.name == "walls" ||
				patch.name == "wall" ||
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
			if (patch.name == "left" ||
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
			else if (patch.name == "right" ||
				patch.name == "outlet")
			{
				bc.emplace_back(make_shared<BoundaryOpen>(patch));
			}
			else if (patch.name == "up" ||
				patch.name == "top" ||
				patch.name == "down" ||
				patch.name == "walls" ||
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


case Blasius:
    {
        double inletRho = phs.cpcv;
        double inletU = 0.1;
        double inletV = 0.0;
        double inletP = 1.0;

        for (const Patch& patch : M.patches)
        {
            if (patch.name == "left" ||
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
            else if (patch.name == "right" ||
                patch.name == "outlet")
            {
                bc.emplace_back(make_shared<BoundaryOpen>(patch));
            }
            else if (patch.name == "up" ||
                patch.name == "top" ||
                patch.name == "top_bottom" ||
                patch.name == "bottom")
            {
                bc.emplace_back(make_shared<BoundarySlip>(patch));
            }
            else if (patch.name == "wall")
            {
                bc.emplace_back(make_shared<BoundaryNonSlip>(patch));
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

case CylinderFlow:
    {
        double inletRho = phs.cpcv;
        double inletU = 0.1;
        double inletV = 0.0;
        double inletP = 1.0;

        for (const Patch& patch : M.patches)
        {
            if (patch.name == "left" ||
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
            else if (patch.name == "right" ||
                patch.name == "outlet")
            {
                bc.emplace_back(make_shared<BoundaryOpen>(patch));
            }
            else if (patch.name == "up" ||
                patch.name == "top" ||
                patch.name == "top_bottom" ||
                patch.name == "bottom")
            {
                bc.emplace_back(make_shared<BoundarySlip>(patch));
            }
            else if (patch.name == "wall")
            {
                bc.emplace_back(make_shared<BoundaryNonSlip>(patch));
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


