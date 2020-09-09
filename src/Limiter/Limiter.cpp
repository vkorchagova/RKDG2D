#include "Limiter.h"
#include <omp.h>

using namespace std;


Limiter::Limiter(
    const Mesh& msh, 
    Solution& sln,
    const Physics& phs,
    const Indicator& ind) 
: mesh(msh), solution(sln), physics(phs), indicator(ind)
{
    newSOL.resize(mesh.nRealCells);
    troubledCells.reserve(mesh.nRealCells);
}


void Limiter::lastHope(std::vector<numvector<double, dimS> >& alpha)
{

/*
    int nCells =  mesh.cells.size();
#pragma omp parallel for schedule (guided) shared(alpha, nCells)
    for (int i = 0; i < nCells; ++i)
	//for (const shared_ptr<Cell>& cell : cells)
    {       
		const shared_ptr<Cell> cell = mesh.cells[i];
		
		for (const shared_ptr<Point>& node : cell->nodes)
        {
            numvector<double, dimPh> res = solution.reconstruct(cell->number, *node);
            double pres = physics.getPressure(res);

            if (res[0] < 0 || res[4] < 0 || pres < 0)
            {

                //cout << "\tnegative values after limitation in cell #" << cell->number << endl;
                //cout << "\t\trho = " << res[RHO] << endl;
                //cout << "\t\trhoU = " << res[RHOU] << endl;
                //cout << "\t\trhoV = " << res[RHOV] << endl;
                //cout << "\t\te = " << res[E] << endl;
                //cout << "\t\tp = " << pres << endl;

                for (int j = 0; j < dimPh; ++j)
                {
                    alpha[cell->number][j*nShapes + 1] = 0.0;
                    alpha[cell->number][j*nShapes + 2] = 0.0;
                }
            }
        }
    }
*/

}


void Limiter::limitSolution()
{        
    //use limiter for all cells
    double t0 = MPI_Wtime();
    troubledCells = indicator.checkDiscontinuities();
    double t1 = MPI_Wtime();
    if (debug) logger << "\t\tLimitSolution.indicator: " << t1 - t0 << endl;
    //cout << "troubledCells.size = " << troubledCells.size() << "\t\t";

    // Limitation stencil
    vector<shared_ptr<Cell>> stencil;
    stencil.reserve(maxPossibleStencilSize);
    newSOL = solution.SOL;

    int nTroubledCells = troubledCells.size();
    int iCell = -1;
    t0 = MPI_Wtime();
//    omp_set_num_threads(NumThreads);
#pragma omp parallel for schedule (guided, maxPossibleStencilSize)\
    shared (nTroubledCells) \
    firstprivate(stencil, iCell) \
    default(none)
    for (int i = 0; i < nTroubledCells; ++i) // here should be range-based cycle but omp doesn't like it :( 
    {
        iCell = troubledCells[i];
		//cout << "iTrCell = " << iCell << endl;
        const shared_ptr<Cell>& cell = mesh.cells[iCell];

        // construct list of cells: cell + neighbours
        stencil = getStencilFor(cell);

        // limit solution
        newSOL[iCell] = limitation(stencil);
    }

    t1 = MPI_Wtime();
    if (debug) logger << "\t\tLimitSolution.forTroubledCells: " << t1 - t0 << endl;

    solution.SOL = newSOL;
    
    t0 = MPI_Wtime();
    lastHope(solution.SOL);
    t1 = MPI_Wtime();
    if (debug) logger << "\t\tLimitSolution.lastHope: " << t1 - t0 << endl;
}
