#include "Solver.h"
#include <iostream>
#include <sstream>
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
    rhs.resize(M.cells.size()); // the same length as SOL
    numFluxes.resize(M.nRealEdges);

    sln.fullSOL.resize(M.nCellsGlob); 
    buf.forFullSOL.resize(M.nCellsGlob * dimS);
    buf.forSolExport.resize(M.nRealCells * dimExp);

	if (myRank == 0)
	{
		buf.forSolExportRecv.resize(M.nCellsGlob * dimExp);
		sln.solToExport.resize(M.nCellsGlob);
	}
        

    ///
    /// preparation for sol collection on proc #0
    ///

    // resize buffers for sol collection on proc #0
    buf.forSendLocalSOL.resize(M.nRealCells * dimS);
    buf.forRecvLocalSOL.resize(M.nRealCells * dimS);

    // resize length/displs for sol collection on proc #0
    buf.nCoeffsPerProc.resize(numProcsTotal);
    buf.mpiDisplSOL.resize(numProcsTotal);

    buf.nSolPerProcExp.resize(numProcsTotal);
    buf.mpiDisplSolExp.resize(numProcsTotal);

    // preparation for gatherv for solution
    if (myRank == 0)
    {
        for (int i = 0; i < numProcsTotal; ++i)
        {
            buf.nCoeffsPerProc[i] = buf.nCellsPerProc[i] * dimS;
            buf.nSolPerProcExp[i] = buf.nCellsPerProc[i] * dimExp;
        }

        buf.mpiDisplSOL[0] = 0;
        buf.mpiDisplSolExp[0] = 0;
        
        for (int i = 1; i < numProcsTotal; ++i)
        {
            buf.mpiDisplSOL[i] = buf.mpiDisplSOL[i-1] + buf.nCoeffsPerProc[i-1];
            buf.mpiDisplSolExp[i] = buf.mpiDisplSolExp[i-1] + buf.nSolPerProcExp[i-1]; 
        }

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
        alpha = B.projection(prb.init, k);
        sln.SOL[k] = alpha;//correctNonOrthoCell(alpha, B.gramian[k]);
    }

} // end setInitialConditions

void Solver::restart(string fileName)
{
    int nCells = M.cells.size();
    sln.SOL.resize(nCells);

    // read full pack of solution
    ifstream reader;
    reader.open(fileName);

    if (!reader.is_open())
    {
        cout << "File " << fileName << " is not found\n";
        exit(0);
    }// if open

    vector<numvector<double, dimS>> fullSOLRestart;
    string line;
    double value;

    while (getline(reader, line))
    {
        int posLine = 0;
        stringstream lineReader(line);
        numvector<double, dimS> currentSOL;

        while (lineReader >> value)
        {
            currentSOL[posLine] = value;
            posLine++;
        }

        fullSOLRestart.push_back(currentSOL);
    }

    reader.close();

    // filter solution values for local MPI proc
    if (numProcsTotal > 1)
    {
        sln.SOL.resize(M.globalCellNumber.size());
        
        for (int iLocal = 0; iLocal < sln.SOL.size(); ++iLocal)
            sln.SOL[iLocal] = fullSOLRestart[M.globalCellNumber[iLocal]];
    }
    else
    {
        sln.SOL = fullSOLRestart;
    }
}


vector<numvector<double, dimS>> Solver::assembleRHS(const std::vector<numvector<double, dimS>>& SOL)
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
        //#pragma omp parallel for \
            shared(myRank, nGP, numFluxes, bcond) \
            private (solLeft, solRight, gpFluxes) \
            default(none)
        //for (const shared_ptr<Edge>& e : bcond->patch.edgeGroup)
        for (int iEdge = 0; iEdge < bcond->patch.edgeGroup.size(); ++iEdge)
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
                solRight = bcond->getSolOuter(solLeft);

                ////if (myRank == 1) cout << "slL = " << solLeft << endl;
                ////if (myRank == 1) cout << "slR = " << solRight << endl;
                
                gpFluxes[iGP] = inverseRotate(flux.evaluate(solLeft, solRight), eNormal);

                ////if (myRank == 1) cout << "; flux: " << gpFluxes[iGP] << endl;
                double nnnn = 0.0;
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
    #pragma omp parallel for  \
         shared(myRank, nGP) \
         firstprivate (solLeft, solRight, gpFluxes) \
         default(none)
    for (int iEdge = M.nEdgesBound; iEdge < M.nRealEdges; ++iEdge)
    {
        ////if (myRank == 1)  cout << "iEdge = " << iEdge <<endl;

        const shared_ptr<Edge>& e = M.edges[iEdge];

        int iCellLeft  = e->neibCells[0]->number;
        int iCellRight = e->neibCells[1]->number;

        ////if (myRank == 1)  cout << iCellLeft << ' ' << iCellRight << endl;


        for (int iGP = 0; iGP < nGP; ++iGP)
        {
            ////if (myRank == 1)  cout << "iGP = " << iGP;// << endl;
            const Point& gPoint = e->gPoints[iGP];
            const Point& eNormal = e->n;
            ////if (myRank == 1)  cout << "gp = " << gPoint << endl;

            solLeft  = rotate(sln.reconstruct(iCellLeft,  gPoint), eNormal);
            solRight = rotate(sln.reconstruct(iCellRight, gPoint), eNormal);

            ////if (myRank == 1) cout << "slL = " << solLeft << endl;
            ////if (myRank == 1) cout << "slR = " << solRight << endl;
            
            gpFluxes[iGP] = inverseRotate(flux.evaluate(solLeft, solRight), eNormal);

            ////if (myRank == 1) cout << "; flux: " << gpFluxes[iGP] << endl;
        }// for GP

        //cout << "-------------" << endl;
        numFluxes[iEdge] = gpFluxes;
        //////if (myRank == 1) cout << "edge #" << iEdge << "; numFlux: " << gpFluxes[0] << ' ' << gpFluxes[1] << endl;

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

numvector<double, dimS> Solver::correctNonOrthoCell(const numvector<double, dimS>& rhs, const vector<vector<double>>& gramian) const
{
    numvector<double, dimS> alphaCorr;

    numvector<double, nShapes> solution(0.0); //for 3 ff!!!

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
#pragma omp parallel /*default(none)*/ \
    shared(alpha, alphaCorr) 
#pragma omp for
    for (int iCell = 0; iCell < M.nRealCells; ++iCell)
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
    ////#pragma omp simd
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
#pragma omp parallel /*default(none)*/ \
    shared(alpha, alphaCorr) 
#pragma omp for
    for (int iCell = 0; iCell < M.nRealCells; ++iCell)
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

	if (numProcsTotal > 1)
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
    if (numProcsTotal > 1)
    {
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
    }
    else
    {
        buf.forFullSOL = buf.forSendLocalSOL;
    }

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


void Solver::collectSolutionForExport()
{
    Point cellCenter;
    for (int iCell = 0; iCell < M.nRealCells; ++iCell)
    {
        cellCenter = M.cells[iCell]->getCellCenter();

        numvector<double, dimPh> solExp = sln.reconstruct(iCell, cellCenter);

        buf.forSolExport[iCell * dimExp] = solExp[0];
        buf.forSolExport[iCell * dimExp + 1] = solExp[1] / solExp[0];
        buf.forSolExport[iCell * dimExp + 2] = solExp[2] / solExp[0];
        buf.forSolExport[iCell * dimExp + 3] = solExp[3] / solExp[0];
        buf.forSolExport[iCell * dimExp + 4] = solExp[4];
        buf.forSolExport[iCell * dimExp + 5] = phs.getPressure(solExp);
    }

	// collect full solution on proc #0
	if (numProcsTotal > 1)
    {
    	MPI_Gatherv(
            &(buf.forSolExport[0]), 
            M.nRealCells * dimExp, 
            MPI_DOUBLE, 
            &(buf.forSolExportRecv[0]), 
            &(buf.nSolPerProcExp[0]), 
            &(buf.mpiDisplSolExp[0]), 
            MPI_DOUBLE, 
            0, 
            MPI_COMM_WORLD
        );
        // sort received solution according to global map

        
    }
    else
    {
        buf.forSolExportRecv = buf.forSolExport;
    }

    for (int i = 0; i < M.nCellsGlob; ++i)
        for (int j = 0; j < dimExp; ++j)
            sln.solToExport[buf.globalMap[i]][j] = buf.forSolExportRecv[i * dimExp + j];
}
