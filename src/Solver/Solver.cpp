#include "Solver.h"
#include <iostream>
#include <omp.h>
#include <mpi.h>

using namespace std;

Solver::Solver( Basis& Bas, Mesh& msh, Solution& soln,
                Problem &prob, Physics& phys, Flux& flx, Buffers& buf) : 
                B(Bas), M(msh), sln(soln), phs(phys), prb(prob), flux(flx), buf(buf)
{

    // // get total number of cells
    // int nCellGlob = 0;
    // MPI_Reduce(&(M.nRealCells), &nCellGlob, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    sln.fullSOL.resize(M.nCellsGlob); 
    buf.forFullSOL.resize(M.nCellsGlob * dimS);

    ///
    /// preparation for sol collection on proc #0
    ///

    // resize buffers for sol collection on proc #0
    buf.forSendLocalSOL.resize(M.nRealCells * dimS);
    buf.forRecvLocalSOL.resize(M.nRealCells * dimS);

    // resize length/displs for sol collection on proc #0
    buf.nCoeffsPerProc.resize(numProcsTotal);
    buf.mpiDisplSOL.resize(numProcsTotal);

    // preparation for gatherv for solution
    if (myRank == 0)
    {
        for (int i = 0; i < numProcsTotal; ++i)
            buf.nCoeffsPerProc[i] = buf.nCellsPerProc[i] * dimS;

        buf.mpiDisplSOL[0] = 0;
        
        for (int i = 1; i < numProcsTotal; ++i)
            buf.mpiDisplSOL[i] = buf.mpiDisplSOL[i-1] + buf.nCoeffsPerProc[i-1]; 

        // cout << "DISPL" << endl;

        // for (int i = 0; i < numProcsTotal; ++i)
        //     cout << buf.mpiDisplSOL[i] << ' ';
        // cout << endl;
    } // if myRank


    ///
    /// preparation for data exchange between boundary procs
    ///

    // resize buffer-sender for proc boundaries
    int nProcCellsTotalInner = 0;
    for (const ProcPatch& p : M.procPatches)
        nProcCellsTotalInner += p.innerCellGroup.size();
    
    buf.forSendBoundSOL.resize(nProcCellsTotalInner * dimS);

    // resize buffer-receiver for proc boundaries
    int nProcCellsTotalOuter = 0;
    for (const ProcPatch& p : M.procPatches)
        nProcCellsTotalOuter += p.cellGroup.size();
    
    buf.forRecvBoundSOL.resize(nProcCellsTotalOuter * dimS);

    // resize len/displs for data exchange
    buf.boundProcSolSizesS.resize(numProcsTotal, 0);
    buf.boundProcSolDisplS.resize(numProcsTotal, 0);

    buf.boundProcSolSizesR.resize(numProcsTotal, 0);
    buf.boundProcSolDisplR.resize(numProcsTotal, 0);

    // prepare len/displs for data send
    for (const ProcPatch& p : M.procPatches)
        for (int i = 0; i < numProcsTotal; ++i) 
            if (i == p.procNum)
            {
                buf.boundProcSolSizesS[i] = p.innerCellGroup.size() * dimS;
                buf.boundProcSolSizesR[i] = p.cellGroup.size() * dimS;
            }

    buf.boundProcSolDisplS[0] = 0;
    buf.boundProcSolDisplR[0] = 0;
    
    for (int i = 1; i < numProcsTotal; ++i)
    {
        buf.boundProcSolDisplS[i] = buf.boundProcSolDisplS[i-1] + buf.boundProcSolSizesS[i-1];
        buf.boundProcSolDisplR[i] = buf.boundProcSolDisplR[i-1] + buf.boundProcSolSizesR[i-1];
    }

    // cout << "rank = " << myRank << "; total mesh size = " << sln.fullSOL.size() << endl;
    // if (myRank == 2)
    // {
    //     cout << "S BOUND SIZES" << endl;

    //     for (int i = 0; i < numProcsTotal; ++i)
    //         cout << buf.boundProcSolSizesS[i] << ' ';
    //     cout << endl;

    //     cout << "S BOUND DISPL" << endl;

    //     for (int i = 0; i < numProcsTotal; ++i)
    //         cout << buf.boundProcSolDisplS[i] << ' ';
    //     cout << endl;

    //     cout << "R BOUND SIZES" << endl;

    //     for (int i = 0; i < numProcsTotal; ++i)
    //         cout << buf.boundProcSolSizesR[i] << ' ';
    //     cout << endl;

    //     cout << "R BOUND DISPL" << endl;

    //     for (int i = 0; i < numProcsTotal; ++i)
    //         cout << buf.boundProcSolDisplR[i] << ' ';
    //     cout << endl;
    // }
    // // get full mapping
}

void Solver::setInitialConditions()
{
    int nCells = M.cells.size();
    numvector<double, dimS> alpha;

    sln.SOL.resize(nCells);

    //#pragma omp parallel for
    for (int k = 0; k < nCells; ++k)
    {
        alpha = projection(prb.init, k);
        sln.SOL[k] = correctNonOrthoCell(alpha, B.gramian[k]);
        //cout << "cell# " << k << "; cfts: " << sln.SOL[k] << endl;
    }

} // end setInitialConditions

void Solver::setDefinedCoefficients(string fileName)
{
    int nCells = M.cells.size();
    sln.SOL.resize(nCells);

    ifstream reader;
    reader.open(fileName);

    if (!reader.is_open())
    {
        cout << "File " << fileName << " is not found\n";
        exit(0);
    }// if open

    numvector<double, dimS> rhs;

    for (int k = 0; k < nCells; ++k)
    {
        for (int j = 0; j < dimS; ++j)
            reader >> sln.SOL[k][j];
    }// for k

    reader.close();
}

numvector<double, dimS> Solver::projection(const std::function<numvector<double, dimPh>(const Point& point)>& foo, int iCell) const
{
    numvector<double, dimS> alpha;    
    numvector<double, dimPh> buffer;
    for (int q = 0; q < nShapes; ++q)
    {
        std::function<numvector<double, dimPh>(const Point&)> f = \
            [&](const Point& p) { return B.phi[q](iCell, p) * foo(p); };

        buffer = integrate(*(M.cells[iCell]), f);
        for (int p = 0; p < dimPh; ++p)
        {
            alpha[p*nShapes + q] = buffer[p];
        }// for p
    }// for shapes

    return alpha;
} // end projection


vector<numvector<double, dimS>> Solver::assembleRHS(const std::vector<numvector<double, dimS>>& SOL)
{
    
    int nCells = M.nRealCells;

    vector<numvector<double, dimS>> rhs(sln.SOL.size()); // the same length as SOL
    vector<vector<numvector<double, dimPh>>> numFluxes(M.nRealEdges);

    double ts0 = 0;
    double te0 = 0;
    
    // 1st step: compute fluxes in gauss points on edges

    int nGP = M.edges[0]->nGP;
    vector<numvector<double, dimPh>> gpFluxes(nGP);
    numvector<double, dimPh> solLeft;
    numvector<double, dimPh> solRight;

    shared_ptr<Edge> e;
    shared_ptr<Cell> cell;
    Point gPoint;
    Point eNormal;

    // for (size_t ind = 0; ind < bc.size(); ++ind)
    //     {
    //         bc[ind]->applyBoundary(sln.SOL);
    //     }

    for (const shared_ptr<Boundary>& bcond : prb.bc)
    {
        for (const shared_ptr<Edge>& e : bcond->patch.edgeGroup)
        {
            int iCellLeft  = e->neibCells[0]->number;

            // cout << iCellLeft << ' ' << iCellRight << endl;
            for (size_t iGP = 0; iGP < nGP; ++iGP)
            {
                // cout << "iGP = " << iGP;// << endl;
                gPoint = e->gPoints[iGP];
                eNormal = e->n;
                // cout << "gp = " << gPoint << endl;

                solLeft  = rotate(sln.reconstruct(iCellLeft,  gPoint), eNormal);
                solRight = bcond->getSolOuter(solLeft);

                // cout << "slL = " << solLeft << endl;
                // cout << "slR = " << solRight << endl;
                
                gpFluxes[iGP] = inverseRotate(flux.evaluate(solLeft, solRight), eNormal);

                // cout << "; flux: " << gpFluxes[iGP] << endl;
            }// for GP

            //cout << "-------------" << endl;
            numFluxes[e->number] = gpFluxes;
            //cout << "edge #" << e->number << "; numFlux: " << gpFluxes[0] << ' ' << gpFluxes[1] << endl;

        }// for bound edges
    } // for bconds 

    for (size_t iEdge = M.nEdgesBound; iEdge < M.nRealEdges; ++iEdge)
    {
        // cout << "iEdge = " << iEdge <<endl;

        e = M.edges[iEdge];

        int iCellLeft  = e->neibCells[0]->number;
        int iCellRight = e->neibCells[1]->number;

        // cout << iCellLeft << ' ' << iCellRight << endl;


        for (size_t iGP = 0; iGP < nGP; ++iGP)
        {
            // cout << "iGP = " << iGP;// << endl;
            gPoint = e->gPoints[iGP];
            eNormal = e->n;
            // cout << "gp = " << gPoint << endl;

            solLeft  = rotate(sln.reconstruct(iCellLeft,  gPoint), eNormal);
            solRight = rotate(sln.reconstruct(iCellRight, gPoint), eNormal);

            // cout << "slL = " << solLeft << endl;
            // cout << "slR = " << solRight << endl;
            
            gpFluxes[iGP] = inverseRotate(flux.evaluate(solLeft, solRight), eNormal);

            // cout << "; flux: " << gpFluxes[iGP] << endl;
        }// for GP

        //cout << "-------------" << endl;
        numFluxes[iEdge] = gpFluxes;
        //cout << "edge #" << iEdge << "; numFlux: " << gpFluxes[0] << ' ' << gpFluxes[1] << endl;

    }// for real edges   


    // 2nd step: compute RHS;

    numvector<double, dimPh> sol;
    numvector<double, dimPh> resV;
    numvector<double, dimS>  res(0.0);
    double gW = 0.0;
    Point nablaPhi;
    //nGP=M.cells[0]->nGP;
    
    for (int iCell = 0; iCell < nCells; ++iCell) // for real (!) cells
    {
        cell = M.cells[iCell];

        // compute internal integral
        res *= 0.0;
        nGP = cell->nGP;

        for (int i = 0; i < nGP; ++i)
        {
            gPoint = cell->gPoints2D[i];
            gW = cell->gWeights2D[i];
            sol = sln.reconstruct(iCell, gPoint);            

            for (int q = 0; q < nShapes; ++q)
            {
                nablaPhi = B.gradPhi[q](iCell, gPoint);
                
                resV = phs.fluxF(sol) * nablaPhi[0] + \
                       phs.fluxG(sol) * nablaPhi[1];

                for (int p = 0; p < dimPh; ++p)
                    res[p * nShapes + q] += resV[p] * gW * cell->J[i];
                 //cout << "cell #" << iCell << "; GP #" << i \
                 //<< "; gradPhi: (" << nablaPhi[0] << ", " << nablaPhi[1] << ")" << endl; 
                 //<< "; res: " << res << endl;     
            }// for shapes


        }// for GP
        rhs[iCell] = res;

        //cout << "cell #" << iCell << "; cellIntegral: " << res << endl;


        // compute boundary integrals

        nGP = M.edges[0]->nGP;
        double sign;
        int iEdge;

        //for (int iEnt = 0; iEnt < cell->nEntities; ++iEnt)
        for (const shared_ptr<Edge>& e : cell->edges)
        {
            res*=0.0;   // ????????????????????

            //iEdge=cell->edges[iEnt]->number;  
            //shared_ptr<Edge> e = M.edges[iEdge];          
            //sign = (iCell == M.edges[iEdge]->neibCells[0]->number) ? 1.0 : -1.0;

            iEdge = e->number;
            sign = (iCell == e->neibCells[0]->number) ? 1.0 : -1.0;

            for (int i = 0; i < nGP; ++i)
            {
                gW = e->gWeights[i]; 
                gPoint = e->gPoints[i]; 

                for (int q = 0; q < nShapes; ++q)
                    for (int p = 0; p < dimPh; ++p)
                        res[p*nShapes + q] += numFluxes[iEdge][i][p] * ( gW * B.phi[q](iCell, gPoint) );
            
            }// for GP

            rhs[iCell] -= res * e->J * sign;
        }// for edges
        
        //cout << "cell #" << iCell << "; RHS[i]: " << rhs[iCell] << endl;


    }// for cells*/ 

    return rhs;
}

numvector<double, dimS> Solver::correctNonOrthoCell(const numvector<double, dimS>& rhs, const vector<vector<double>>& gramian) const
{
    numvector<double, dimS> alphaCorr;

    vector<double> solution(nShapes); //for 3 ff!!!

    //cout << "cell# " << iCell << "; cfts: " << alpha << "; G: " << B.gramian[iCell] << endl;
    for (int iSol = 0; iSol < dimPh; ++iSol)
    {
        // solve slae
        solution[0] = rhs[iSol*nShapes];


        if (nShapes == 3)
        {
            //cout << "\tiSol = " << iSol << endl;

            solution[2] = (rhs[iSol*nShapes + 2] * gramian[0][0] - rhs[iSol*nShapes + 1] * gramian[1][0]) \
                / (gramian[1][1] * gramian[0][0] - gramian[1][0] * gramian[1][0]);
            solution[1] = (rhs[iSol*nShapes + 1] - solution[2] * gramian[1][0]) \
                  / (gramian[0][0]);
        }

        //set solution to appropriate positions
        for (int i = 0; i < nShapes; ++i)
        {
            //cout << "i = " << i << "; iAlpha = " << i + iSol*nShapes << endl; 
            alphaCorr[i + iSol*nShapes] = solution[i];
        }

        //cout << "the next sol..." << endl;
    }

    //cout << "result in fun = " << alphaCorr << endl;
    //cout << "the next cell..." << endl;

    return alphaCorr;
}


vector<numvector<double, dimS>> Solver::correctNonOrtho(const vector<numvector<double, dimS>>& alpha) const
{
    vector<numvector<double, dimS>> alphaCorr(alpha.size());

    //cout << "size of alpha " << alpha.size() << endl;
    //cout << "size of alphaCorr " << alphaCorr.size() << endl;
    for (size_t iCell = 0; iCell < M.nRealCells; ++iCell)
    {
        //cout << "iCell in common  = " << iCell << endl;
        alphaCorr[iCell] = correctNonOrthoCell(alpha[iCell], B.gramian[iCell]);
        //cout << "result = " << alphaCorr[iCell] << endl;
    }// for cells

    return alphaCorr;
}


numvector<double, dimS> Solver::correctPrevIterCell(const numvector<double, dimS>& alphaCorr, const vector<vector<double>>& gramian) const
{
    numvector<double, dimS> alpha(0.0);
    for (int iSol = 0; iSol < dimPh; ++iSol)
    {
        alpha[iSol*nShapes] += alphaCorr[iSol*nShapes];

        for (int i = 1; i < nShapes; ++i)
        {
    //#pragma omp simd
            for (int j = 1; j <=i; ++j)
                alpha[i + iSol*nShapes] += gramian[i-1][j-1] * alphaCorr[iSol*nShapes + j];

            for (int j = i + 1; j < nShapes; ++j)
                alpha[i + iSol*nShapes] += gramian[j-1][i-1] * alphaCorr[iSol*nShapes + j];

        }
    }// for variables

    return alpha;
}


vector<numvector<double, dimS>> Solver::correctPrevIter(const vector<numvector<double, dimS>>& alpha) const
{
    vector<numvector<double, dimS>> alphaCorr(alpha.size());

    //cout << "size of alpha " << alpha.size() << endl;
    //cout << "size of alphaCorr " << alphaCorr.size() << endl;
    for (size_t iCell = 0; iCell < M.nRealCells; ++iCell)
    {
        //cout << "iCell in common  = " << iCell << endl;
        alphaCorr[iCell] = correctPrevIterCell(alpha[iCell], B.gramian[iCell]);
        //cout << "result = " << alphaCorr[iCell] << endl;
    }// for cells

    return alphaCorr;
}


void Solver::dataExchange()
{
    // 1st step: fill bound send buffers and resize recv buffers

    for (const ProcPatch& p : M.procPatches)
    {
        int curDispl = buf.boundProcSolDisplS[p.procNum];
        // if (myRank == 2)
        //     cout << "proc num = " << p.procNum << endl;

        for (int i = 0; i < p.innerCellGroup.size(); ++i)
        {
            int iCell = p.innerCellGroup[i]->number;

            for (int j = 0; j < dimS; ++j)
                buf.forSendBoundSOL[curDispl + i*dimS + j] = sln.SOL[iCell][j];
            
            // if (myRank == 2)
            // {
            //     cout << "\nbufferize sol on cell #" << M.globalCellNumber[iCell] << endl;
            //     for (int j = 0; j < dimS; ++j)
            //         cout << sln.SOL[iCell][j] << ' ';
            //     cout << endl;
            // }
        }
    }

    // if (myRank == 2)
    // {
    //     cout << "buf.forSendBoundSOL" << endl;
    //     int stopper = 0;
    //     for (int i = 0; i < buf.forSendBoundSOL.size(); ++i)
    //     {    
    //         cout << buf.forSendBoundSOL[i] << ' ';
    //         if (stopper == dimS -1)
    //         {
    //             cout << endl;
    //             stopper = 0;
    //         }
    //         else
    //             ++stopper;
    //     }
    // }

    // 2nd step: scatter data for all procs

    MPI_Alltoallv(
        &(buf.forSendBoundSOL[0]),      /*who send*/
        &(buf.boundProcSolSizesS[0]),   /*how much*/
        &(buf.boundProcSolDisplS[0]),   /*displs*/
        MPI_DOUBLE,
        &(buf.forRecvBoundSOL[0]),      /*who recv*/
        &(buf.boundProcSolSizesR[0]),   /*how much*/
        &(buf.boundProcSolDisplR[0]),   /*displs*/
        MPI_DOUBLE,
        MPI_COMM_WORLD
    );

    //MPI_Barrier(MPI_COMM_WORLD);

    // 3rd step: sort data after receiving

    for (const ProcPatch& p : M.procPatches)
    {
        int curDispl = buf.boundProcSolDisplR[p.procNum];

        // if (myRank == 2)
        //     cout << "proc num = " << p.procNum << endl;

        for (int i = 0; i < p.cellGroup.size(); ++i)
        {
            int iCell = p.cellGroup[i]->number;
            
            for (int j = 0; j < dimS; ++j)
                sln.SOL[iCell][j] = buf.forRecvBoundSOL[curDispl + i*dimS + j];
        }
    }



    // if (myRank == 2)
    // {
    //     cout << "buf.forRecvBoundSOL" << endl;
    //     int stopper = 0;
    //     for (int i = 0; i < buf.forRecvBoundSOL.size(); ++i)
    //     {    
    //         cout << buf.forRecvBoundSOL[i] << ' ';
    //         if (stopper == dimS - 1)
    //         {
    //             cout << endl;
    //             stopper = 0;
    //         }
    //         else
    //             ++stopper;
    //     }
    // }

    // if (myRank == 2)
    // {
    //     cout << "sol " << M.globalCellNumber[4] << " on proc 2\n" << sln.SOL[4] << endl;
    // }

    // if (myRank == 3)
    // {
    //     cout << "sol "<< M.globalCellNumber[2] << " on proc 3\n" << sln.SOL[2] << endl;
    // }
}

void Solver::collectSolution()
{
    // put solution on each proc into buffer

    for (int i = 0; i < M.nRealCells; ++i)
        for (int j = 0; j < dimS; ++j)
            buf.forSendLocalSOL[i*dimS + j] = sln.SOL[i][j];

    // if (myRank == 0)
    // {
    //     cout << "sendBuf for local sol on proc 0" << endl;
    //     for (int i = 0; i < M.nRealCells; ++i)
    //     {   
    //         for (int j = 0; j < dimS; ++j)
    //             cout << buf.forSendLocalSOL[i*dimS + j] << ' ';
    //         cout << endl;
    //     }
    // }



    // collect full solution on proc #0
    MPI_Gatherv(
        &(buf.forSendLocalSOL[0]), 
        M.nRealCells*dimS, 
        MPI_DOUBLE, 
        &(buf.forFullSOL[0]), 
        &(buf.nCoeffsPerProc[0]), 
        &(buf.mpiDisplSOL[0]), 
        MPI_DOUBLE, 
        0, 
        MPI_COMM_WORLD
    );

    //MPI_Barrier(MPI_COMM_WORLD);

    // sort received solution according to global map

    for (int i = 0; i < M.nCellsGlob; ++i)
        for (int j = 0; j < dimS; ++j)
            sln.fullSOL[buf.globalMap[i]][j] = buf.forFullSOL[i*dimS + j];

    // if (myRank == 0)
    // {
    //     cout << "buf.forFullSOL for local sol on proc 0" << endl;
    //     for (int i = 0; i < M.nCellsGlob; ++i)
    //     {   
    //         for (int j = 0; j < dimS; ++j)
    //             cout << buf.forFullSOL[i*dimS + j] << ' ';
    //         cout << endl;
    //     }

    //     // cout << "\nfullSOL for local sol on proc 0" << endl;
    //     // for (int i = 0; i < M.nCellsGlob; ++i)
    //     // {   
    //     //     for (int j = 0; j < dimS; ++j)
    //     //         cout << sln.fullSOL[i][j] << ' ';
    //     //     cout << endl;
    //     // }
    // }
    
}