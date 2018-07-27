#include "Problem.h"
#include <iostream>

#include "Patch.h"
#include "EdgeBoundary.h"

using namespace std;

// ------------------ Constructors & Destructors ----------------

//Problem::Problem(const std::vector<numvector<double, 5 * nShapes> > &al) : alpha(al)
Problem::Problem(string caseName, const Time& t) : time(t)
{
    setInitialConditions(caseName);
} // end constructor by mesh

Problem::~Problem()
{

}


// ------------------ Private class methods --------------------

// ------------------ Public class methods --------------------



void Problem::setInitialConditions(string caseName)
{
    // Function for initial condition
    //double rho0 = 0.5;

    function<double(const Point& r)> initRho;
    function<double(const Point& r)> initP;
    function<double(const Point& r)> initU;
    function<double(const Point& r)> initV;

    if (caseName == "SodX")
    {
        cpcv = 1.4;

        initRho = [](const Point& r) { return (r.x() < 0.0) ? 1.0 : 0.125; };
        initP   = [](const Point& r) { return (r.x() < 0.0) ? 1.0 : 0.1;  };
        initU   = [](const Point& r) { return (r.x() < 0.0) ? 0.0 : 0.0; };
        initV   = [](const Point& r) { return 0.0; };
    }
    else if (caseName == "SodY")
    {
        cpcv = 1.4;

        initRho = [](const Point& r) { return (r.x() < 0) ? 1.0 : 0.125; };
        initP   = [](const Point& r) { return (r.x() < 0) ? 1.0 : 0.1;  };
        initV   = [](const Point& r) { return (r.x() < 0) ? 0.75 : 0.0; };
        initU   = [](const Point& r) { return 0.0; };
    }
    else if (caseName == "SodDiag")
    {
        cpcv = 1.4;

        initRho = [](const Point& r) { return ((r.x() + r.y()) < 1.0) ? 1.0 : 0.125; };
        initP   = [](const Point& r) { return ((r.x() + r.y()) < 1.0) ? 1.0 : 0.1;  };
        initV   = [](const Point& r) { return 0.0; };
        initU   = [](const Point& r) { return 0.0; };
    }
    else if (caseName == "SodCircle")
    {
        cpcv = 1.4;

        initRho = [](const Point& r) { return (r.length() <= 0.4) ? 1.0 : 0.125; };
        initP   = [](const Point& r) { return (r.length() <= 0.4) ? 1.0 : 0.1;  };
        initV   = [](const Point& r) { return 0.0; };
        initU   = [](const Point& r) { return 0.0; };
    }
    else if (caseName == "Woodward")
    {
        cpcv = 1.4;

        initRho = [](const Point& r) { return 1.0; };
        initP   = [](const Point& r) { return (r.x() <= -0.4) ? 1000.0 : ((r.x() >= 0.4) ? 100.0 : 0.01);  };
        initV   = [](const Point& r) { return 0.0; };
        initU   = [](const Point& r) { return 0.0; };
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
    else if (caseName == "acousticPulse")
    {
        cpcv = 1.4;

        initRho = [](const Point& r) { return 1.0 + 1e-6*exp( - 40.0*sqr(r.x() )- 40.0*sqr(r.y() )); };
        initP   = [&](const Point& r) { return initRho(r) / cpcv;  };
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
        return numvector<double, 5> { initRho(r), initU(r), initV(r), 0.0, initP(r) / (cpcv - 1.0) + 0.5 * initRho(r) * (sqr(initU(r)) + sqr(initV(r)))};
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
        caseName == "Noh" || \
        caseName == "123" || \
        caseName == "acousticPulse")
    {
        shared_ptr<BoundaryOpen> bOpen = make_shared<BoundaryOpen>();
        shared_ptr<BoundarySlip> bSlip = make_shared<BoundarySlip>();

        //bc = {bOpen, bOpen};
        bc = {bSlip, bSlip};
        //bc = {bSlip, bSlip, bSlip, bSlip};
        //bc = {bOpen, bOpen, bOpen, bOpen};
    }
    else if (caseName == "forwardStep")
    {
        shared_ptr<BoundaryOpen> bOpen = make_shared<BoundaryOpen>();
        shared_ptr<BoundarySlip> bSlip = make_shared<BoundarySlip>();
        shared_ptr<BoundaryConstant> bConst = \
                make_shared<BoundaryConstant>(numvector<double,5>({1.0, -3.0, 0.0, 0.0, 6.286}));

        bc = {bConst, bOpen, bSlip, bSlip};
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

//// RKDG methods

void Problem::setAlpha(const std::vector<numvector<double, 5 * nShapes> >& a)
{
    alpha = a;
}

double Problem::getPressure(const numvector<double, 5>& sol) const
{
    // uncomment for LEE
   numvector<double,5> initfun = init(Point({0.0,0.0}));

   double rho0 = initfun[0];
   double p0 = initfun[4] * (cpcv - 1);

   return p0 * pow(sol[0] / rho0 , cpcv);

    // end uncomment for LEE

    double magRhoU2 = sqr(sol[1]) + sqr(sol[2]) + sqr(sol[3]);

    double p = (cpcv - 1.0)*(sol[4] - 0.5*magRhoU2 / sol[0]);

    return p;
} // end getPressure

double Problem::c(const numvector<double, 5>& sol) const
{
    double c2 = cpcv * getPressure(sol) / sol[0];

//    if (c2 < 0)
//        cout << "!!!!" << endl;

    return sqrt( c2 );
} // end c for cell


double Problem::c_av(const numvector<double, 5>& solOne, const numvector<double, 5>& solTwo) const
{
    double semiRho = 0.5*(solOne[0] + solTwo[0]);
    double semiP =   0.5*(getPressure(solOne) + getPressure(solTwo));

    return sqrt( cpcv * semiP / semiRho);
} // end c for edge

numvector<double, 5> Problem::lambdaF_Roe(const numvector<double, 5>& solOne, const numvector<double, 5>& solTwo) const
{
    double sqrtRhoLeft = sqrt(solOne[0]);
    double sqrtRhoRight = sqrt(solTwo[0]);



    double sumSqrtRho = sqrtRhoLeft + sqrtRhoRight;

    double u_av = ( solOne[1] / sqrtRhoLeft + solTwo[1] / sqrtRhoRight ) / sumSqrtRho;
    double h_av = ( (solOne[4] + getPressure(solOne)) / sqrtRhoLeft + (solTwo[4] + getPressure(solTwo)) / sqrtRhoRight ) / sumSqrtRho;

    double c_av = sqrt( (cpcv - 1) * (h_av - 0.5 * sqr(u_av)) );

//    if (getPressure(solOne) < 0 || getPressure(solTwo) < 0)
//        cout << "kkk";

    return {u_av - c_av, u_av, u_av, u_av, u_av + c_av};

}

numvector<double, 5> Problem::lambdaF_Einfeldt(const numvector<double, 5>& solOne, const numvector<double, 5>& solTwo) const
{
    double cLeft  = c(solOne);
    double cRight = c(solTwo);

    double uLeft = solOne[1] / solOne[0];
    double uRight = solTwo[1] / solTwo[0];

    double sqrtRhoLeft = sqrt(solOne[0]);
    double sqrtRhoRight = sqrt(solTwo[0]);
    double sumSqrtRho = sqrtRhoLeft + sqrtRhoRight;

    double eta2 = 0.5 * sqrtRhoLeft * sqrtRhoRight / sqr(sqrtRhoLeft + sqrtRhoRight);
    double u_av =  ( solOne[1] / sqrtRhoLeft + solTwo[1] / sqrtRhoRight ) / sumSqrtRho;
    double c_av = sqrt( (sqrtRhoLeft * sqr(cLeft) + sqrtRhoRight * sqr(cRight)) / sumSqrtRho +\
                        eta2 * sqr(uRight - uLeft));

    return {u_av - c_av, u_av, u_av, u_av, u_av + c_av};

}

numvector<double, 5> Problem::lambdaF_Toro(const numvector<double, 5>& solOne, const numvector<double, 5>& solTwo) const
{
    double uLeft = solOne[1] / solOne[0];
    double uRight = solTwo[1] / solTwo[0];

    double cLeft  = c(solOne);
    double cRight = c(solTwo);

    double pLeft = getPressure(solOne);
    double pRight = getPressure(solTwo);

    double pvrs = 0.5 * (pLeft + pRight) + \
            0.125 * (uLeft - uRight) * (solOne[0] + solTwo[0]) * (cLeft + cRight);

    double pStar = max(0.0, pvrs);

    double qLeft = pStar > pLeft ? \
                    sqrt( 1.0 + 0.5 * (cpcv + 1.0) * (pStar / pLeft - 1.0) / cpcv) : \
                    1.0;

    double qRight = pStar > pRight ? \
                    sqrt( 1.0 + 0.5 * (cpcv + 1.0) * (pStar / pRight - 1.0) / cpcv) : \
                    1.0;

    return {uLeft - qLeft * cLeft, uLeft, uLeft, uLeft, uRight + qRight * cRight};
}

numvector<double, 5> Problem::lambdaF_semisum(const numvector<double, 5>& solOne, const numvector<double, 5>& solTwo) const
{
    double u = 0.5*(solOne[1] / solOne[0] + solTwo[1] / solTwo[0]);
    double soundSpeed = c_av(solOne, solTwo);

    return { u - soundSpeed, u, u, u, u + soundSpeed};
} // end lambdaF

numvector<double, 5> Problem::lambdaF(const numvector<double, 5>& solOne, const numvector<double, 5>& solTwo) const
{
    return lambdaF_Roe(solOne,solTwo);
} // end lambdaF



numvector<double, 5> Problem::fluxF(const numvector<double, 5>& sol) const
{
    double u = sol[1] / sol[0];
    double p = getPressure(sol);

    return { sol[1], u*sol[1] + p, u*sol[2], u*sol[3], (sol[4] + p)*u };
} // end fluxF

numvector<double, 5> Problem::fluxG(const numvector<double, 5>& sol) const
{
    double v = sol[2] / sol[0];
    double p = getPressure(sol);

    return { sol[2], v*sol[1], v*sol[2] + p, v*sol[3], (sol[4] + p)*v };
} // end fluxG


pair<numvector<numvector<double, 5>, 5>, numvector<numvector<double, 5>, 5>> Problem::getL(const numvector<double, 5>& sol) const
{
    double cS = c(sol);
    double rho = sol[0];
    double u = sol[1] / rho;
    double v = sol[2] / rho;

    double B1 = (cpcv - 1.0) / sqr(cS);
    double B2 = 0.5 * B1 * (sqr(u) + sqr(v));

    pair<numvector<numvector<double, 5>, 5>, numvector<numvector<double, 5>, 5>> res;

    res.first =
    {
        // Lx
        { 0.5 * (B2 + u/cS), -0.5 * (B1 * u + 1.0/cS), -0.5 * B1 * v,  0.0, 0.5 * B1},
        {                -v,                      0.0,           1.0,  0.0,      0.0},
        {               0.0,                      0.0,           0.0,  1.0,      0.0},
        {          1.0 - B2,                   B1 * u,        B1 * v,  0.0,      -B1},
        { 0.5 * (B2 - u/cS), -0.5 * (B1 * u - 1.0/cS), -0.5 * B1 * v,  0.0, 0.5 * B1}
    };

    res.second =
    {
        // Ly
        { 0.5 * (B2 + v/cS), -0.5 * (B1 * u), -0.5 * (B1 * v + 1.0/cS),  0.0, 0.5 * B1},
        {                 u,            -1.0,                      0.0,  0.0,      0.0},
        {               0.0,             0.0,                      0.0,  1.0,      0.0},
        {          1.0 - B2,          B1 * u,                   B1 * v,  0.0,      -B1},
        { 0.5 * (B2 - v/cS), -0.5 * (B1 * u), -0.5 * (B1 * v - 1.0/cS),  0.0, 0.5 * B1}
    };

    return res;
}


pair<numvector<numvector<double, 5>, 5>, numvector<numvector<double, 5>, 5>> Problem::getR(const numvector<double, 5>& sol) const
{
    double cS = c(sol);
    double rho = sol[0];
    double u = sol[1] / rho;
    double v = sol[2] / rho;
    double p = getPressure(sol);

    double H = (sol[4] + p) / rho;
    double magU2 = 0.5 * (sqr(u) + sqr(v));

    pair<numvector<numvector<double, 5>, 5>, numvector<numvector<double, 5>, 5>> res;

    res.first =
    {
        // Rx
        {        1.0,  0.0,  0.0,   1.0,        1.0 },
        {     u - cS,  0.0,  0.0,     u,     u + cS },
        {          v,  1.0,  0.0,     v,          v },
        {        0.0,  0.0,  1.0,   0.0,        0.0 },
        { H - cS * u,    v,  0.0, magU2, H + cS * u }
    };

    res.second =
    {
        // Ry
        {        1.0,  0.0,  0.0,   1.0,        1.0 },
        {          u, -1.0,  0.0,     u,          u },
        {     v - cS,  0.0,  0.0,     v,     v + cS },
        {        0.0,  0.0,  1.0,   0.0,        0.0 },
        { H - cS * v,   -u,  0.0, magU2, H + cS * v },
    };

    return res;
}

numvector<numvector<double, 5>, 5> Problem::getL(const numvector<double, 5>& sol, const Point& n) const
{
    double cS = c(sol);
    double rho = sol[0];
    double u = sol[1] / rho;
    double v = sol[2] / rho;

    double B1 = (cpcv - 1.0) / sqr(cS);
    double B2 = 0.5 * B1 * (sqr(u) + sqr(v));

    numvector<numvector<double, 5>, 5> res;

    res =
    {
        { 0.5 * (B2 + (u*n.x() + v*n.y())/cS), -0.5 * (B1 * u + n.x()/cS), -0.5 * (B1 * v + n.y()/cS),  0.0, 0.5 * B1},
        {                   u*n.y() - v*n.x(),                     -n.y(),                      n.x(),  0.0,      0.0},
        {                            1.0 - B2,                     B1 * u,                     B1 * v,  0.0,      -B1},
        {                                 0.0,                        0.0,                        0.0,  1.0,      0.0},
        { 0.5 * (B2 - (u*n.x() + v*n.y())/cS), -0.5 * (B1 * u - n.x()/cS), -0.5 * (B1 * v - n.y()/cS),  0.0, 0.5 * B1}
    };

    return res;
}


numvector<numvector<double, 5>, 5> Problem::getR(const numvector<double, 5>& sol, const Point& n) const
{
    double cS = c(sol);
    double rho = sol[0];
    double u = sol[1] / rho;
    double v = sol[2] / rho;
    double p = getPressure(sol);

    double H = (sol[4] + p) / rho;
    double magU2 = 0.5 * (sqr(u) + sqr(v));

    numvector<numvector<double, 5>, 5> res;

    res =
    {
        {                          1.0,                0.0,   1.0,  0.0,                          1.0 },
        {               u - cS * n.x(),             -n.y(),     u,  0.0,                 u + cS*n.x() },
        {               v - cS * n.y(),              n.x(),     v,  0.0,                 v + cS*n.y() },
        {                          0.0,                0.0,   0.0,  1.0,                          0.0 },
        { H - cS * (u*n.x() + v*n.y()), -u*n.y() + v*n.x(), magU2,  0.0, H + cS * (u*n.x() + v*n.y()) }
    };

    return res;
}









