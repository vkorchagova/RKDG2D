#include "SolverHybrid.h"
#include <iostream>
#include <sstream>
#include <omp.h>
#include <mpi.h>

using namespace std;

SolverHybrid::SolverHybrid( Basis& Bas, Mesh& msh, Solution& soln,
                Problem &prob, Physics& phys, Flux& flx1, Flux& flx2, Buffers& buf) : 
                Solver(Bas, msh, soln, prob, phys, flx1, buf), flux1(flx1), flux2(flx2)
{
}


void SolverHybrid::getHybridFlux( 
    const numvector<double, dimPh>& solLeft, 
    const numvector<double, dimPh>& solRight,
    const Point& eNormal,
    numvector<double, dimPh>& flux
)
{
    double momMag = sqrt((solLeft[1] - solRight[1])*(solLeft[1] - solRight[1]) + (solLeft[2] - solRight[2])*(solLeft[2] - solRight[2]));

    Point n1 = momMag > 1e-3 ? Point({(solLeft[1] - solRight[1]) / momMag, (solLeft[2] - solRight[2]) / momMag }) : eNormal;
    Point n2 ({n1.y(), -n1.x()}); // which direction for n2?

    double alpha1 = eNormal * n1;
    double alpha2 = eNormal * n2;

    if (alpha1 < 0.0)
    {
        alpha1 *= -1;
        n1 *= -1;
    }

    if (alpha2 < 0.0)
    {
        alpha2 *= -1;
        n2 *= -1;
    }

    if (alpha1 < 1e-6)
    {
        alpha1 = 0.0;
        alpha2 = 1.0;
    }

    if (alpha2 < 1e-6)
    {
        alpha2 = 0.0;
        alpha1 = 1.0;
    }

    numvector<double, dimPh> sol1n1 = rotate(solLeft, n1);
    numvector<double, dimPh> sol2n1 = rotate(solRight, n1);

    numvector<double, dimPh> sol1n2 = rotate(solLeft, n2);
    numvector<double, dimPh> sol2n2 = rotate(solRight, n2);

    flux = alpha1 * inverseRotate(flux1.evaluate(sol1n1, sol2n1), n1) + alpha2 * inverseRotate(flux2.evaluate(sol1n2, sol2n2), n2);
    // cout << "alpha1 = " << alpha1 << "; alpha2 = " << alpha2 << "; flux = " << flux << endl;
}


vector<numvector<double, dimS>> SolverHybrid::assembleRHS(const std::vector<numvector<double, dimS>>& SOL)
{   
    int nCells = M.nRealCells;

    double t0, t1;
    
    // 1st step: compute fluxes in gauss points on edges

    int nGP = M.edges[0]->nGP;
    vector<numvector<double, dimPh>> gpFluxes(nGP);
    numvector<double, dimPh> solLeft;
    numvector<double, dimPh> solRight;

    // cout << phs.c(sln.reconstruct(M.cells[0]->number, M.cells[0]->getCellCenter()));
    // exit(1);

    ///--------------------------------------------------------------------------------
    t0 = MPI_Wtime();

    for (const shared_ptr<Boundary>& bcond : prb.bc)
    {
        int nEdgesPatch = bcond->patch.edgeGroup.size();
        #pragma omp parallel for \
            shared(nGP, numFluxes, bcond, nEdgesPatch) \
            firstprivate (solLeft, solRight, gpFluxes) \
            default(none)
        //for (const shared_ptr<Edge>& e : bcond->patch.edgeGroup)
        for (int iEdge = 0; iEdge < nEdgesPatch; ++iEdge)
        {
            const shared_ptr<Edge>& e = bcond->patch.edgeGroup[iEdge];
            ////if (myRank == 1) cout << e->number << endl;
            int iCellLeft  = e->neibCells[0]->number;
            ////if (myRank == 1) cout << "-------------" << endl;

            ////if (myRank == 1) cout << iCellLeft << ' '  << endl;
            for (size_t iGP = 0; iGP < nGP; ++iGP)
            {
               ////if (myRank == 1) cout << "iGP = " << iGP;// << endl;
                Point& gPoint = e->gPoints[iGP];
                Point& eNormal = e->n;
                ////if (myRank == 1) cout << "gp = " << gPoint << endl;

                solLeft  = rotate(sln.reconstruct(iCellLeft,  gPoint), eNormal);
                solRight = bcond->getSolOuter(solLeft, eNormal);

                ////if (myRank == 1) cout << "slL = " << solLeft << endl;
                ////if (myRank == 1) cout << "slR = " << solRight << endl;
                
                getHybridFlux(inverseRotate(solLeft, eNormal), inverseRotate(solRight, eNormal), eNormal, gpFluxes[iGP]);

                ////if (myRank == 1) cout << "; flux: " << gpFluxes[iGP] << endl;
            }// for GP

            
            numFluxes[e->number] = gpFluxes;
            ////if (myRank == 1) cout << "edge #" << e->number << "; numFlux: " << gpFluxes[0] << ' ' << gpFluxes[1] << endl;

        }// for bound edges
    } // for bconds 
    t1 = MPI_Wtime();
    if (debug) logger << "\t\teBound.numfluxes: " << t1 - t0 << endl;
    ////if (myRank == 1) cout << "end bound edges" << endl;


    ///--------------------------------------------------------------------------------
    t0 = MPI_Wtime();

    //omp_set_num_threads(NumThreads);
    #pragma omp parallel for schedule (guided)  \
         shared(myRank, nGP, cout) \
         firstprivate (solLeft, solRight, gpFluxes) \
         default(none)
    for (int iEdge = M.nEdgesBound; iEdge < M.nRealEdges; ++iEdge)
    {
        ////if (myRank == 1)  cout << "iEdge = " << iEdge <<endl;

        const shared_ptr<Edge>& e = M.edges[iEdge];



        int iCellLeft  = e->neibCells[0]->number;
        // cout << iCellLeft << endl;
        // cout << e->nodes[0]->x() << ' ' 
        //      << e->nodes[0]->y() << ' '
        //      << e->nodes[1]->x() << ' '
        //      << e->nodes[1]->y() << endl;
        int iCellRight = e->neibCells[1]->number;

        


        for (int iGP = 0; iGP < nGP; ++iGP)
        {
            ////if (myRank == 1)  cout << "iGP = " << iGP;// << endl;
            const Point& gPoint = e->gPoints[iGP];
            const Point& eNormal = e->n;
            ////if (myRank == 1)  cout << "gp = " << gPoint << endl;

            solLeft  = sln.reconstruct(iCellLeft,  gPoint);
            solRight = sln.reconstruct(iCellRight, gPoint);

            ////if (myRank == 1) cout << "slL = " << solLeft << endl;
            ////if (myRank == 1) cout << "slR = " << solRight << endl;
            
            getHybridFlux(solLeft, solRight, eNormal, gpFluxes[iGP]);

            ////if (myRank == 1) cout << "; flux: " << gpFluxes[iGP] << endl;
            //if ( (iCellLeft == 50) || (iCellRight == 50) || (iCellLeft == 1562) || (iCellRight == 1562))
            //{
            //    cout << "-------" << endl;
            //    cout << "iCellLeft = " << iCellLeft << "; iCellRight = " << iCellRight << endl;
            //    cout << "solLeft = " << solLeft << endl;
            //    cout << "solRight = " << solRight << endl;
            //    cout << "numFlux: " << gpFluxes[iGP] << endl;
            //}
        }// for GP

        //cout << "-------------" << endl;
        numFluxes[iEdge] = gpFluxes;
        //////if (myRank == 1) 

        
    }// for real edges   
    t1 = MPI_Wtime();
    if (debug) logger << "\t\teInner.numfluxes: " << t1 - t0 << endl;
   
    ///--------------------------------------------------------------------------------
    // 2nd step: compute RHS;

    t0 = MPI_Wtime();

    #pragma omp parallel default(none) \
     shared(myRank, nCells) \
     private(nGP)
    {
    numvector<double, dimPh> sol;
    numvector<double, dimPh> resV;
    numvector<double, dimS>  res(0.0);
    double gW = 0.0;
    Point nablaPhi;

    #pragma omp for
    for (int iCell = 0; iCell < nCells; ++iCell) // for real (!) cells
    {
        const shared_ptr<Cell>& cell = M.cells[iCell];
        res *= 0.0;

        ///----------------------------------------------------
        // compute internal integral

        nGP = cell->nGP;
        double coef = 0.0; // aux var for optimization

        for (int i = 0; i < nGP; ++i)
        {
            const Point& gPoint = cell->gPoints2D[i];
            gW = cell->gWeights2D[i];
            sol = sln.reconstruct(iCell, gPoint);            
            coef = gW * cell->J[i];

            for (int q = 0; q < nShapes; ++q)
            {
                nablaPhi = B.gradPhi[q](iCell, gPoint);
                
                resV = phs.fluxF(sol) * nablaPhi[0] + \
                       phs.fluxG(sol) * nablaPhi[1];

                resV += prb.source(sol, gPoint) * B.phi[q](iCell, gPoint);

                for (int p = 0; p < dimPh; ++p)
                    res[p * nShapes + q] += resV[p] * coef;    
            }// for shapes


        }// for GP
        rhs[iCell] = res;

        //cout << "cell #" << iCell << "; cellIntegral: " << res << endl;


        ///----------------------------------------------------
        // compute boundary integrals
        nGP = M.edges[0]->nGP;

        double sign;
        int iEdge;

        for (const shared_ptr<Edge>& e : cell->edges)
        {
            res *= 0.0;   // ZEROFICATION!

            iEdge = e->number;
            sign = (iCell == e->neibCells[0]->number) ? 1.0 : -1.0;
            coef = 0.0;

            for (int i = 0; i < nGP; ++i)
            {
                gW = e->gWeights[i]; 
                const Point& gPoint = e->gPoints[i]; 

                for (int q = 0; q < nShapes; ++q)
                {
                    coef = gW * B.phi[q](iCell, gPoint);
                    for (int p = 0; p < dimPh; ++p)
                        res[p*nShapes + q] += numFluxes[iEdge][i][p] * coef;
                }// for shapes
            }// for GP

            rhs[iCell] -= res * e->J * sign;

        }// for edges
        
        //cout << "cell #" << iCell << "; RHS[i]: " << rhs[iCell] << endl;

    }// for cells*/ 
    }// omp parallel
    t1 = MPI_Wtime();
    if (debug) logger << "\t\trhs.compute: " << t1 - t0 << endl;

    return rhs;
}