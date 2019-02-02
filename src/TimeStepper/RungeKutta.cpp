#include "RungeKutta.h"

using namespace std;

RungeKutta::RungeKutta(int o,  Basis& b, Solver& s, Solution& ss, std::vector<std::shared_ptr<Boundary>>& bond, Limiter& l, TimeControl& t) : TimeStepper(o, b, s, ss, bond, l, t)
{
    /// Initialization of the parameters
	// Number of stages a.k.a. length of all the arrays
	nStages = order; // OK for RK order <= 4

	// The stages array sizing
    k.resize(nStages);
	//Arr.resize(nStages);

	// Cfts arrays sizing
    alpha.resize(nStages);
    beta.resize(nStages);
    for (int i = 0; i < nStages; ++i)
        beta[i].resize(nStages);

	// Initialization of the RK-cfts
    setButcherTable();

}

void RungeKutta::setButcherTable()
{
    switch (order)
    {
        case 1:

            beta[0][0] = 1.0;
            break;

        case 2:

            alpha[0] = 0.0;
            alpha[1] = 1.0;

            //beta[0][0] = 0.5;
            //beta[1][0] = 0.0;
            //beta[1][1] = 1.0;

            beta[0][0] = 1.0;
            beta[1][0] = 0.5;
            beta[1][1] = 0.5;

            break;

        case 3:

            alpha[1] = 0.5;
            alpha[2] = 1.0;

            beta[0][0] = 0.5;
            beta[1][1] = 1.0;

            beta[2][0] = 0.1666666666666667;
            beta[2][1] = 0.6666666666666667;
            beta[2][2] = 0.1666666666666667;

            break;

        case 4:

            alpha[1] = 0.5;
            alpha[2] = 0.5;
            alpha[3] = 1.0;

            beta[0][0] = 0.5;
            //beta[1][0] = 0.0;
            beta[1][1] = 0.5;
            //beta[2][0] = 0.0;
            //beta[2][1] = 0.0;
            beta[2][2] = 1.0;

            beta[3][0] = 0.1666666666666667;
            beta[3][1] = 0.3333333333333333;
            beta[3][2] = 0.3333333333333333;
            beta[3][3] = 0.1666666666666667;

            break;

        default:
            cout << "Runge --- Kutta method of order " \
                 << order \
                 << "is not implemented. " \
                 << "Please change order to 1, 2, 3 or 4." \
                 << endl;
            exit(1);

    }
}

void RungeKutta::Tstep()
{
	/// Final preparations
    double t = T.getTime();
	double tau = T.getTau();

							///!!! SOL_aux must be equal SOL at this point !!!
    vector<numvector<double, dimS>> lhs    = sln.SOL;
    vector<numvector<double, dimS>> lhsOld = slv.correctPrevIter(sln.SOL);

    //for(int iCell=0; iCell<sln.SOL.size(); ++iCell)\
    //    cout << "cell#" << iCell << "; SOL: " << sln.SOL[iCell] << endl;
    //cout << endl;

    //cout << "OK" << endl;
	/// The very step of the RK method
    for (int i = 0; i < nStages; ++i)
    {
        T.updateTime(t + alpha[i]*tau);   // ??? Is it necessary?

        // MPI exchange between neib procs

        slv.dataExchange();
        
         //if (myRank == 0)
         //{
         //   cout << "sln sol before rhs" << endl;
         ///   for (int p = 0; p < sln.SOL.size(); ++p)
         //        cout << p << ' ' << sln.SOL[p] << endl;
         //}

        // assemble rhs of SODE
        k[i] = slv.assembleRHS(sln.SOL);

        // if (myRank == 0)
        // {
        //    cout << "k" << i << endl;
        //    for (int p = 0; p < k[i].size(); ++p)
        //         cout << p << ' ' << k[i][p] << endl;
        // }


        // RK step
        lhs = lhsOld;

        for (int j = 0; j <= i; ++j)
           lhs = lhs + k[j] * beta[i][j] * tau; 

        // if (myRank == 0)
        // {
        //    cout << "lhs after rhs" << endl;
        //    for (int p = 0; p < lhs.size(); ++p)
        //         cout << p << ' ' << lhs[p] << endl;
        // }

        // remember about non-ortho basis!
        sln.SOL = slv.correctNonOrtho(lhs);
	
	    // if (myRank == 0)
        // {
        //    cout << "sln sol after non ortho" << endl;
        //    for (int p = 0; p < sln.SOL.size(); ++p)
        //         cout << p << ' ' << sln.SOL[p] << endl;
        // }

        // for(int iCell = 0; iCell < sln.SOL.size(); ++iCell)\
        //    cout << "cell#" << iCell << "; SOL: " << sln.SOL[iCell] << endl;
        // cout << endl;
		
		// limit solution
        slv.dataExchange();
        lmt.limit(sln.SOL);

    }// for stages
    
	/// Updating the solution and time step
	T.updateTimeStep(slv.MaxSpeed); // !!! IT MUST BE DEFINED ALREADY !!!
	
    //for(int iCell=0; iCell<16; ++iCell)\
        cout << "cell#" << iCell << "; SOL: " << sln.SOL[iCell] << endl;
    //cout << endl;

    //time.updateTime(tOld);

    return;
}