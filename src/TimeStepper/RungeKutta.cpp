#include "RungeKutta.h"
 
using namespace std;

RungeKutta::RungeKutta(int o,  Basis& b, Solver& s, Solution& ss, Limiter& l, TimeControl& t) : TimeStepper(o, b, s, ss, l, t)
{
    // Initialization of the parameters
    nStages = order; // OK for RK order <= 4
    
    if (order < 0) // just for modifications
        nStages = 3;

    if (order < -2) // just for new modifications
        nStages = 4;

    // Cfts arrays sizing 
    alpha.resize(nStages, 0.0);
    beta.resize(nStages);
    for (int i = 0; i < nStages; ++i)
        beta[i].resize(nStages, 0.0);

    // Initialization of the RK-cfts
    setButcherTable();
}


void RungeKutta::setButcherTable()
{
    switch (order)
    {
        case 1:

            alpha[0] = 1.0;
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
        
            alpha[2] = 1.0;

            beta[0][0] = 1.0;
            beta[1][0] = 0.25;
            beta[1][1] = 0.25;

            beta[2][0] = 0.1666666666666667;
            beta[2][1] = 0.1666666666666667;
            beta[2][2] = 0.6666666666666667;
            
            break;

        case -1:
	    
	    alpha[2] = 1.0;

            beta[0][0] = 0.5;
            
            beta[1][0] = 0.5;
            beta[1][1] = 0.5;

            beta[2][0] = 0.3333333333333333;
            beta[2][1] = 0.3333333333333333;
            beta[2][2] = 0.3333333333333333;
            
            break;
            
        case -2:

            alpha[2] = 1.0;

            beta[0][0] = 0.3333333333333333;
            
            beta[1][0] = 0.3333333333333333;
            beta[1][1] = 0.3333333333333333;

            beta[2][0] = 0.25;
            beta[2][1] = 0.0;
            beta[2][2] = 0.75;

            break;
        
        case -3:

            alpha[3] = 1.0;

            beta[0][0] = 0.3333333333333333;
            
            beta[1][0] = 0.3333333333333333;
            beta[1][1] = 0.3333333333333333;

            beta[2][0] = 0.3333333333333333;
            beta[2][1] = 0.3333333333333333;
            beta[2][2] = 0.3333333333333333;
            
            beta[3][0] = 0.25;
            beta[3][1] = 0.25;
            beta[3][2] = 0.25;
            beta[3][3] = 0.25;

            break;
            
        case -4:
        {
            alpha[3] = 1.0;
            
            double pp = 0.6;
            double aa = 6.0/11.0;
            double bb = 56.0/111.0;
            double gg = 31.0/216.0;

            beta[0][0] = pp;
            
            beta[1][0] = pp * (1.0 - aa);
            beta[1][1] = pp * (1.0 - aa);

            beta[2][0] = pp * (1.0 - aa) * (1.0 - bb);
            beta[2][1] = pp * (1.0 - aa) * (1.0 - bb);
            beta[2][2] = pp * (1.0 - bb);

            beta[3][0] = pp * (1.0 - aa) * (1.0 - bb) * (1.0 - gg);
            beta[3][1] = pp * (1.0 - aa) * (1.0 - bb) * (1.0 - gg);
            beta[3][2] = pp * (1.0 - bb) * (1.0 - gg);
            beta[3][3] = pp * (1.0 - gg);

            break;
        }

        case -5:
        {
            alpha[3] = 1.0;
            
            double pp = 0.528;//0.57132;
            double aa = 0.2514970059880235;//0.4611604448659367;
            double bb = 0.601766535035007;//0.5334643755038199;
            double gg = 0.05036722250048653;//0.11119560206959966;

            beta[0][0] = pp;
            
            beta[1][0] = pp * (1.0 - aa);
            beta[1][1] = pp * (1.0 - aa);

            beta[2][0] = pp * (1.0 - aa) * (1.0 - bb);
            beta[2][1] = pp * (1.0 - aa) * (1.0 - bb);
            beta[2][2] = pp * (1.0 - bb);

            beta[3][0] = pp * (1.0 - aa) * (1.0 - bb) * (1.0 - gg);
            beta[3][1] = pp * (1.0 - aa) * (1.0 - bb) * (1.0 - gg);
            beta[3][2] = pp * (1.0 - bb) * (1.0 - gg);
            beta[3][3] = pp * (1.0 - gg);

            break;
        }

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
                 << "Use order -1 for RK2BC1 or -2 for RK2BC2." \
                 << "Use order -3 for RK4ORD2 or -4 for RK4ORD3." \
                 << endl;
            exit(1);
    }
}


void RungeKutta::Tstep()
{
    double t0, t1;

    /// Final preparations
    double t = T.getTime();
    double tau = T.getTau();

    /// rhs for RK studies
    std::vector<std::vector<numvector<double, dimS>>> k;
    // The stages array sizing
    k.resize(nStages);

    vector<numvector<double, dimS>> lhs    = sln.SOL;

    t0 = MPI_Wtime();
    vector<numvector<double, dimS>> lhsOld = slv.correctPrevIter(sln.SOL);
    t1 = MPI_Wtime();
    if (debug) logger << "\tslv.correctPrevIter(): " << t1 - t0 << endl;

    //for(int iCell=0; iCell<sln.SOL.size(); ++iCell)\
    //    cout << "cell#" << iCell << "; SOL: " << sln.SOL[iCell] << endl;
    //cout << endl;

    //cout << "OK" << endl;
    /// The very step of the RK method
    for (int iStage = 0; iStage < nStages; ++iStage)
    {
        T.updateTimeValueRK(alpha[iStage]*tau);   // ??? Is it necessary?

        // MPI exchange between neib procs
        t0 = MPI_Wtime();
        slv.dataExchange();
        t1 = MPI_Wtime();
        if (debug) logger << "\tslv.dataExchange(): " << t1 - t0 << endl;

        //cout << myRank << "__after data exchange" << endl;
        
         // if (myRank == 0)
         // {
         //   cout << "sln sol before rhs" << endl;
         // /   for (int p = 0; p < sln.SOL.size(); ++p)
         //        cout << p << ' ' << sln.SOL[p] << endl;
         // }

        // assemble rhs 
        t0 = MPI_Wtime();
        k[iStage] = slv.assembleRHS(sln.SOL);
        t1 = MPI_Wtime();
        if (debug) logger << "\tslv.assembleRHS(): " << t1 - t0 << endl;


        // if (myRank == 0)
        // {
        //    cout << "k" << i << endl;
        //    for (int p = 0; p < k[i].size(); ++p)
        //         cout << p << ' ' << k[i][p] << endl;
        // }


        // RK step
        t0 = MPI_Wtime();


        #pragma omp parallel for \
            shared(lhs, lhsOld) \
            default(none)
        for (int iCell = 0; iCell < lhs.size(); ++iCell)
            for (int iSol = 0; iSol < dimS; ++iSol)
                lhs[iCell][iSol] = lhsOld[iCell][iSol];

        
        for (int jStage = 0; jStage <= iStage; ++jStage)
        {
            double rkCoeff = beta[iStage][jStage] * tau;

            #pragma omp parallel for \
                shared(lhs, k, rkCoeff, iStage, jStage) \
                default(none) 
            for (int iCell = 0; iCell < lhs.size(); ++iCell)
                for (int iSol = 0; iSol < dimS; ++iSol)
                    lhs[iCell][iSol] += k[jStage][iCell][iSol] * rkCoeff; 
        }

        //lhs = lhsOld;

        //for (int j = 0; j <= iStage; ++j)
        //    lhs += k[j] * (beta[iStage][j] * tau);

        t1 = MPI_Wtime();
        if (debug) logger << "\tupdateRKstep: " << t1 - t0 << endl;

        // if (myRank == 0)
        // {
        //    cout << "lhs after rhs" << endl;
        //    for (int p = 0; p < lhs.size(); ++p)
        //         cout << p << ' ' << lhs[p] << endl;
        // }

        // remember about non-ortho basis!
        t0 = MPI_Wtime();
        sln.SOL = slv.correctNonOrtho(lhs);
        t1 = MPI_Wtime();
        if (debug) logger << "\tslv.correctNonOrtho(): " << t1 - t0 << endl;
    
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
        t0 = MPI_Wtime();
        slv.dataExchange();
        t1 = MPI_Wtime();
        if (debug) logger << "\tslv.dataExchange(): " << t1 - t0 << endl;

        t0 = MPI_Wtime();
        lmt.limitSolution();
        t1 = MPI_Wtime();
        if (debug) logger << "\tslmt.limit(): " << t1 - t0 << "\n\t-----" << endl;

        //cout << "\t END RK STAGE #" << iStage << endl;

    }// for stages

    return;
}