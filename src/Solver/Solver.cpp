#include "Solver.h"
#include <iostream>
#include <omp.h>

using namespace std;

Solver::Solver( Basis& Bas, Mesh& msh, Solution& soln,
                Problem &prob, Physics& phys, Flux& flx) : 
                B(Bas), M(msh), sln(soln), phs(phys), prb(prob), flux(flx)
{
    
}

void Solver::setInitialConditions()
{
    int nCells = M.cells.size();
    numvector<double, dimS> alpha;

    sln.SOL.resize(nCells);
    sln.SOL_aux.resize(nCells);

    //#pragma omp parallel for
    for (int k = 0; k < nCells; ++k)
    {
        alpha = projection(prb.init, k);
        sln.SOL[k] = correctNonOrthoCell(alpha, B.gramian[k]);
        //cout << "cell# " << k << "; cfts: " << sln.SOL[k] << endl;
    }
    
    cout << "OK" << endl;

} // end setInitialConditions

void Solver::setDefinedCoefficients(string fileName)
{
    ifstream reader;
    reader.open(fileName);

    if (!reader.is_open())
    {
        cout << "File " << fileName << " is not found\n";
        exit(0);
    }// if open

    int nCells = M.cells.size();

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

    for (size_t iEdge = 0; iEdge < M.nRealEdges; ++iEdge)
    {
        //cout << "iEdge = " << iEdge <<endl;

        shared_ptr<Edge> e = M.edges[iEdge];

        int iCellLeft  = e->neibCells[0]->number;
        int iCellRight = e->neibCells[1]->number;

        //cout << iCellLeft << ' ' << iCellRight << endl;


        for (size_t iGP = 0; iGP < nGP; ++iGP)
        {
            //cout << "iGP = " << iGP;// << endl;

            solLeft  = sln.reconstruct(iCellLeft,  e->gPoints[iGP]);
            solRight = sln.reconstruct(iCellRight, e->gPoints[iGP]);

            //cout << "slL = " << solLeft << endl;
            //cout << "slR = " << solRight << endl;
            
            gpFluxes[iGP] = flux.evaluate(solLeft, solRight, e->n);

            //cout << "; flux: " << gpFluxes[iGP] << endl;
        }// for GP

        //cout << "-------------" << endl;
        numFluxes[iEdge] = gpFluxes;
        //cout << "edge #" << iEdge << "; numFlux: " << res << endl;

    }// for edges   



    // 2nd step: compute RHS;

    numvector<double, dimPh> sol;
    numvector<double, dimPh> resV;
    numvector<double, dimS> res(0.0);
    double gW = 0.0;
    Point pt;
    Point nablaPhi;
    //nGP=M.cells[0]->nGP;
    
    for (int iCell = 0; iCell < nCells; ++iCell) // for real (!) cells
    {
        shared_ptr<Cell> cell = M.cells[iCell];

        // compute internal integral
        res *= 0.0;
        nGP = cell->nGP;

        for (int i = 0; i < nGP; ++i)
        {
            pt = cell->gPoints2D[i];
            gW = cell->gWeights2D[i];
            sol = sln.reconstruct(iCell, pt);            

            for (int q = 0; q < nShapes; ++q)
            {
                nablaPhi = B.gradPhi[q](iCell, pt);
                
                resV = phs.fluxF(sol) * nablaPhi[0] + \
                       phs.fluxG(sol) * nablaPhi[1];

                for (int p = 0; p < dimPh; ++p)
                    res[p * nShapes + q] += resV[p] * gW * cell->J[i]; // CHECK!!!
            //cout << "cell #" << iCell << "; GP #" << i \
                 //<< "; gradPhi: (" << nablaPhi[0] << ", " << nablaPhi[1] << ")" << endl; 
                 //<< "; res: " << res << endl;     
            }// for shapes


        }// for GP
        rhs[iCell] = res;

        if (iCell > M.nRealCells)
            cout << "WARNING: iCell = " << iCell << "more than mesh size!!" << endl;   

        //cout << "cell #" << iCell << "; cellIntegral: " << res << endl;


        // compute boundary integrals

        nGP = M.edges[0]->nGP;
        double sign;
        int iEdge;

        for (int iEnt = 0; iEnt < cell->nEntities; ++iEnt)
        {
            res*=0.0;   // ????????????????????

            iEdge=cell->edges[iEnt]->number;  
            shared_ptr<Edge> e = M.edges[iEdge];          
            sign = (iCell == M.edges[iEdge]->neibCells[0]->number) ? 1.0 : -1.0;

            for (int i = 0; i < nGP; ++i)
            {
                gW = e->gWeights[i]; 
                pt = e->gPoints[i]; 

                for (int q = 0; q < nShapes; ++q)
                    for (int p = 0; p < dimPh; ++p)
                        res[p*nShapes + q] += numFluxes[iEdge][i][p] * ( gW * B.phi[q](iCell, pt) );
            
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