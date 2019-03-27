#include "LimiterRiemannWENOS.h"

using namespace std;

const int maxPossibleStencilSize = dimPh;
const double g = 0.001;

#pragma omp threadprivate(maxPossibleStencilSize, g)


numvector<double, dimS> limitP
(
    const vector<shared_ptr<Cell>>& stenc, 
    const vector<numvector<double, dimS>>& pInv,
    const Solution& solution,
    vector<double> gamma,
    vector<numvector<double, dimPh>> beta,
    vector<numvector<double, dimPh>> w,
    vector<numvector<double, dimPh>> wTilde,
    numvector<double, dimPh> wSum,
    vector<numvector<double, dimPh>> uMean
)
{   // just for limiting of px and py

    int n = pInv.size(); 
    int nCells = stenc.size(); 

    

    /// p polynomas
    std::vector<numvector<double, dimS>> p = pInv;


    // get mean values of linear functions

    uMean.resize(nCells);

    for (size_t k = 0; k < nCells; ++k)
        //uMean[k] = sln.reconstructSolution(cells[0]->getCellCenter());
        uMean[k] = solution.reconstruct(stenc[k]->number, stenc[0]->getCellCenter(), p[k]);

    // get coeffs for polynoms p = a0 + a1 * (x - xc) + a2 * (y - yc)

    for (size_t k = 0; k < nCells; ++k)
        for (int i = 0; i < dimPh; ++i)
            for (int j = 0; j < nShapes; ++j)
                p[k][i*nShapes + j] *= solution.B.phiCoeffs[k][j];
/*
cout << "after phiCoeffs" << endl;
cout << "cell #" << cell->number << ":\n";
    
for (size_t k = 0; k < nCells; ++k) \
            cout << "p #" << k << ": " << p[k] << endl;
*/

    for (size_t k = 0; k < nCells; ++k)
        for (int i = 0; i < dimPh; ++i)
            p[k][i*nShapes] += - uMean[k][i] + uMean[0][i];

    // get linear weights

    gamma.resize(nCells);

    gamma[0] = 1.0 - (nCells - 1) * g;

    for (size_t k = 1; k < nCells; ++k)
        gamma[k] = g;

    // get smoothness indicators and nonlinear weights

    beta.resize(nCells);
    wTilde.resize(nCells);
    w.resize(nCells);

    for (size_t k = 0; k < nCells; ++k)
        for (int j = 0; j < dimPh; ++j)
        {
            beta[k][j] = stenc[0]->getArea() * stenc[k]->getArea() *
                (sqr(p[k][j*nShapes + 1]) + sqr(p[k][j*nShapes + 2]));
            wTilde[k][j] = gamma[k] * (1.0 / sqr(beta[k][j] + 1e-6));
        }



    //wSum = {0.0, 0.0, 0.0, 0.0, 0.0};
    for (size_t k = 0; k < dimPh; ++k)
        wSum[k] = 0.0;


    for (int j = 0; j < dimPh; ++j)
        for (size_t k = 0; k < nCells; ++k)
            wSum[j] += wTilde[k][j];

//        //cout << wSum << endl;

    for (size_t k = 0; k < nCells; ++k)
        for (int j = 0; j < dimPh; ++j)
            w[k][j] = wTilde[k][j] / wSum[j];

//        cout << "----\n num tr cell = " << iCell << endl;

//        cout << "beta:\n";
//        for (size_t k = 0; k < nCells; ++k)
//        {
//            cout << "cell no = " << k << endl;
//            cout << beta[k] << endl;
//        }

//        cout << "wtilde:\n";
//        for (size_t k = 0; k < nCells; ++k)
//        {
//            cout << "cell no = " << k << endl;
//            cout << wTilde[k] << endl;
//        }

//        cout << "w:\n";
//        for (size_t k = 0; k < nCells; ++k)
//        {
//            cout << "cell no = " << k << endl;
//            cout << w[k] << endl;
//        }

    numvector<double, dimS> res (0.0);

    for (int i = 0; i < dimPh; ++i)
    {
        res[i*nShapes] =  p[0][i*nShapes];

        for (size_t k = 0; k < nCells; ++k)
        {
            //res[i*nShapes]     += w[k][i] * p[k][i*nShapes];
            res[i*nShapes + 1] += w[k][i] * p[k][i*nShapes + 1];
            res[i*nShapes + 2] += w[k][i] * p[k][i*nShapes + 2];
        }
    }

    return res;
}

void LimiterRiemannWENOS::limit(vector<numvector<double, dimS>>& alpha)
{
    int n = alpha.size();

    double t0, t1;
    t0 = MPI_Wtime();

    //troubledCells = indicator.checkDiscontinuities();
    troubledCells.resize(n);
#pragma omp parallel for \
 shared(troubledCells, n) \
 default(none)
    for (int i = 0; i < n; ++i)
        troubledCells[i] = i; 

    /// Linear weights
    std::vector<double> gamma;

    /// smoothness indicators
    std::vector<numvector<double, dimPh>> beta;

    /// nonlinear weights
    std::vector<numvector<double, dimPh>> w;
    std::vector<numvector<double, dimPh>> wTilde;
    numvector<double, dimPh> wSum;

    /// mean values
    std::vector<numvector<double, dimPh>> uMean;

    gamma.reserve(maxPossibleStencilSize);
    beta.reserve(maxPossibleStencilSize);
    w.reserve(maxPossibleStencilSize);
    wTilde.reserve(maxPossibleStencilSize);
    uMean.reserve(maxPossibleStencilSize);

    // p polynoms

    vector<numvector<double, dimS>> p;
    vector<numvector<double, dimS>> pInv;
    vector<numvector<double, dimS>> pNew;

    p.reserve(maxPossibleStencilSize);
    pInv.reserve(maxPossibleStencilSize);
    pNew.reserve(maxPossibleStencilSize);

    // stencil
    vector<shared_ptr<Cell>> stenc;
    stenc.reserve(maxPossibleStencilSize); 

    // limit solution in troubled cells

    numvector<numvector<double, dimPh>, dimPh>  L;
    numvector<numvector<double, dimPh>, dimPh>  R;
    numvector <double, dimS> pLim;
    numvector<double, dimPh> solMean;

#pragma omp parallel for \
shared(alpha, cout) \
private(uMean, beta, gamma, w, wTilde, wSum, p, pInv, pNew, stenc, L, R, pLim, solMean) \
default(none)
    for (int iTroubled = 0; iTroubled < troubledCells.size(); ++iTroubled)
    //for (int iCell : troubledCells)
    {
        int iCell = troubledCells[iTroubled];

        shared_ptr<Cell> cell = cells[iCell];

        // construct list of cells: cell + neighbours

        stenc = { cell };
        stenc.insert(stenc.end(), cell->neibCells.begin(), cell->neibCells.end());

        int nCells = stenc.size();

        // get p polynoms for each cell

        p.resize(nCells);
        pInv.resize(nCells);

        for (int k = 0; k < nCells; ++k)
            p[k] = alpha[stenc[k]->number];

         //cout << "before sol mean" << endl;
            
        solMean = solution.reconstruct(cell->number,cell->getCellCenter());

        //cout << "after sol mean" << endl;
        
        // get pNew polynoms using each edge normal
        for (const shared_ptr<Edge> e : cell->edges)
        {            
            //cout << "processing edge #" << e->number << endl;
            if (e->neibCells.size() == 2)
            {
                // correct normal direction: outside related to troubled cell
                Point n = (( *e->nodes[0] - cell->getCellCenter()) * e->n > 0.0) ? e->n : Point(-e->n);

                L = physics.getL(solMean, n);
                R = physics.getR(solMean, n);

                // get Riemann invariants along this direction
                for (size_t k = 0; k < stenc.size(); ++k)
                    pInv[k] = conservativeToRiemann(p[k], L);


                // limit Riemann invariants
                pLim = limitP(stenc, pInv, solution, gamma, beta, w ,wTilde, wSum, uMean);


                // project limited solution to conservative variables
                pNew.push_back( riemannToConservative(pLim, R) );
            }
        }

        double sumArea = 0.0;

        for (int i = 1; i < nCells; ++i)
            sumArea += stenc[i]->getArea();

        // construct limited solution
        numvector<double, dimS> res (0.0);
        for (int i = 0; i < nCells-1; ++i)
            res += pNew[i] * (stenc[i+1]->getArea() / sumArea);

        



        function<numvector<double, dimPh>(const Point& r)> foo = [=](const Point& r) \
        {
            numvector<double, dimPh> sum (0.0);

            for (int i = 0; i < dimPh; ++i)
                sum[i] += (res[i*nShapes] + \
                           res[i*nShapes + 1] * (r.x() - cell->getCellCenter().x()) + \
                           res[i*nShapes + 2] * (r.y() - cell->getCellCenter().y()));


            return sum;
        };

         alphaNew[iCell] = solution.B.projection(foo, cell->number);
    }

    alpha = alphaNew;
    
    t1 = MPI_Wtime();
    if (debug) logger << "\t\tLimiterRiemannWENOS.forTroubledCells: " << t1 - t0 << endl;

    alpha = alphaNew;

    t0 = MPI_Wtime();
    //lastHope(alpha);
    t1 = MPI_Wtime();
    if (debug) logger << "\t\tLimiterRiemannWENOS.lastHope: " << t1 - t0 << endl;
}



numvector<double, dimS> LimiterRiemannWENOS::conservativeToRiemann
(
    const numvector<double, dimS>& alpha, 
    const numvector<numvector<double, dimPh>, dimPh>& L
) const
{
    numvector<double, dimS> res (0.0);

    for (int iSol = 0; iSol < dimPh; ++iSol )
        for (int jSol = 0; jSol < dimPh; ++jSol )
            for (int j = 0; j < nShapes; ++j )
                res[iSol * nShapes + j] +=  L[iSol][jSol] * alpha[jSol * nShapes + j];

    return res;
}

numvector<double, dimS> LimiterRiemannWENOS::riemannToConservative
(
    const numvector<double, dimS>& alpha, 
    const numvector<numvector<double, dimPh>, dimPh>& R
) const
{
    numvector<double, dimS> res (0);

    for (int iSol = 0; iSol < dimPh; ++iSol )
        for (int jSol = 0; jSol < dimPh; ++jSol )
            for (int j = 0; j < nShapes; ++j )
                res[iSol * nShapes + j] += R[iSol][jSol] * alpha[jSol * nShapes + j];

    return res;

}



