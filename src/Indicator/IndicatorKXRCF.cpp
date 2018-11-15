#include "IndicatorKXRCF.h"
#include <omp.h>

using namespace std;

bool debugFlux = false;

const double threshold = 1e-12;

#pragma omp threadprivate(debugFlux)

numvector<double, 4> IndicatorKXRCF::massFlux(const Edge& edge, const Cell& cell) const
{
    double f[2];

    Point n = Point({0.0, 0.0});

    Point n0 = *edge.nodes[0];
    Point n1 = *edge.nodes[1];

    numvector<double, 5>  solOne, solTwo;


        
    try // find internal edges
	{
        // for internal edges
        const EdgeInternal& edgeInt = dynamic_cast<const EdgeInternal&>(edge);

        // neighbour cell could not be equal to considered cell
        Cell& neib = (&cell == edge.neibCells[0].get()) ? *edge.neibCells[1] : *edge.neibCells[0];

        // correct normal position: outside related to considered cell
        n = ((n0 - cell.getCellCenter()) * edge.n > 0.0) ? edge.n : Point(-edge.n);

        // get normal velocity component in nodes of edge... averaged by Roe!
        solOne = rotate(cell.reconstructSolution(n0), n);
        solTwo = rotate(neib.reconstructSolution(n0), n);
        f[0] = solOne[1] / sqrt(solOne[0]) + solTwo[1] / sqrt(solTwo[0]);

        solOne = rotate(cell.reconstructSolution(n1), n);
        solTwo = rotate(neib.reconstructSolution(n1), n);
        f[1] = solOne[1] / sqrt(solOne[0]) + solTwo[1] / sqrt(solTwo[0]);

        // 4 variants of flux integral (only through the part of edge with inward flux)
        // we have linear functions -> use trapezia formulae for integration


        // no entrance
        if ((f[0] >= -threshold) && (f[1] >= -threshold))
        {
            if (debugFlux)
            {
//#pragma omp master
                cout << "no entrance\n";
            }
            return { 0.0, 0.0, 0.0, 0.0 };
        }
        double h = 0.0;
        double intRho = 0.0;
        double intE = 0.0;
        double intP = 0.0;

        //entrance at n[0]
        if (f[0] < -threshold)
		{
            Point nZero = (f[1] <= -threshold) ? \
                           n1 : \
                           Point({ (f[1] * n0.x() - f[0] * n1.x())/ (f[1] - f[0]), (f[1] * n0.y() - f[0] * n1.y())/ (f[1] - f[0])});

            h = (nZero - n0).length();

            intRho = 0.5*h*(
                        cell.reconstructSolution(n0, 0) + cell.reconstructSolution(nZero, 0) -
                        neib.reconstructSolution(n0, 0) - neib.reconstructSolution(nZero, 0)
                  );

            intE = 0.5*h*(
                        cell.reconstructSolution(n0, 4) + cell.reconstructSolution(nZero, 4) -
                        neib.reconstructSolution(n0, 4) - neib.reconstructSolution(nZero, 4)
                  );
            intP = 0.5*h*(
                        problem.getPressure(cell.reconstructSolution(n0)) + problem.getPressure(cell.reconstructSolution(nZero)) -
                        problem.getPressure(neib.reconstructSolution(n0)) - problem.getPressure(neib.reconstructSolution(nZero))
                  );

            if (debugFlux)
            {
//#pragma omp master
                cout << "\n----\nu0 < 0 int \n -- \n";
                cout << "h = " << h << endl;
                cout << "nZero = " << nZero << endl;
                cout << "mySol zero rho = "   << cell.reconstructSolution(nZero, 0) << endl;
                cout << "neibSol zero rho = " << neib.reconstructSolution(nZero, 0) << endl;
                cout << "mySol n1 rho = "   << cell.reconstructSolution(n1, 0) << endl;
                cout << "neibSol n1 rho = " << neib.reconstructSolution(n1, 0) << endl;
                cout << "mySol zero e = "     << cell.reconstructSolution(nZero, 4) << endl;
                cout << "neibSol zero e = "   << neib.reconstructSolution(nZero, 4) << endl;
                cout << "mySol n0 rho = "     << cell.reconstructSolution(n0, 0) << endl;
                cout << "neibSol n0 rho = "   << neib.reconstructSolution(n0, 0) << endl;
                cout << "mySol n0 e = "       << cell.reconstructSolution(n0, 4) << endl;
                cout << "neibSol n0 e = "     << neib.reconstructSolution(n0, 4) << endl;

                cout << "int rho = "          << intRho << endl;
                cout << "int e = "            << intE << endl;
            }
		}
        else // no entrance at n0, entrance at n1
		{
            Point nZero = Point({ (f[1] * n0.x() - f[0] * n1.x())/ (f[1] - f[0]), (f[1] * n0.y() - f[0] * n1.y())/ (f[1] - f[0])});

            h = (n1 - nZero).length();

            intRho = 0.5*h*(
                        cell.reconstructSolution(n1, 0) + cell.reconstructSolution(nZero, 0) -
                        neib.reconstructSolution(n1, 0) - neib.reconstructSolution(nZero, 0)
                  );

            intE = 0.5*h*(
                        cell.reconstructSolution(n1, 4) + cell.reconstructSolution(nZero, 4) -
                        neib.reconstructSolution(n1, 4) - neib.reconstructSolution(nZero, 4)
                  );

            intP = 0.5*h*(
                        problem.getPressure(cell.reconstructSolution(n1)) + problem.getPressure(cell.reconstructSolution(nZero)) -
                        problem.getPressure(neib.reconstructSolution(n1)) - problem.getPressure(neib.reconstructSolution(nZero))
                  );

            if (debugFlux)
            {
//#pragma omp master
                cout << "\n----\nu0 > 0, u1 < 0 int \n -- \n";
                cout << "h = " << h << endl;
                cout << "nZero = " << nZero << endl;
                cout << "mySol zero rho = "   << cell.reconstructSolution(nZero, 0) << endl;
                cout << "neibSol zero rho = " << neib.reconstructSolution(nZero, 0) << endl;
                cout << "mySol zero e = "     << cell.reconstructSolution(nZero, 4) << endl;
                cout << "neibSol zero e = "   << neib.reconstructSolution(nZero, 4) << endl;
                cout << "mySol n1 rho = "     << cell.reconstructSolution(n1, 0) << endl;
                cout << "neibSol n1 rho = "   << neib.reconstructSolution(n1, 0) << endl;
                cout << "mySol n1 e = "       << cell.reconstructSolution(n1, 4) << endl;
                cout << "neibSol n1 e = "     << neib.reconstructSolution(n1, 4) << endl;

                cout << "int rho = "          << intRho << endl;
                cout << "int e = "            << intE << endl;
            }
		}

        return
        {
            intRho,
            intE,
            intP,
            h
        };

    } // end try
    catch (...) // other edges are boundaries
	{
        // for edge on boundary
        const EdgeBoundary& edgeBound = dynamic_cast<const EdgeBoundary&>(edge);

        n = edge.n;

        // get normal velocity component in nodes of edge ... averaged by Roe!
        solOne = rotate(cell.reconstructSolution(n0), n);
        solTwo = rotate(edgeBound.bc->applyBoundary(cell.reconstructSolution(n0)), n);
        f[0] = solOne[1] / sqrt(solOne[0]) + solTwo[1] / sqrt(solTwo[0]);

        solOne = rotate(cell.reconstructSolution(n1), n);
        solTwo = rotate(edgeBound.bc->applyBoundary(cell.reconstructSolution(n1)), n);
        f[1] = solOne[1] / sqrt(solOne[0]) + solTwo[1] / sqrt(solTwo[0]);




        if ((f[0] >= -threshold) && (f[1] >= -threshold))
        {
            if (debugFlux)
            {
//#pragma omp master
                cout << "no entrance\n";
            }
            return { 0.0, 0.0, 0.0, 0.0 };
        }

        double h = 0.0;
        double intRho = 0.0;
        double intE = 0.0;
        double intP = 0.0;

        if (f[0] < -threshold)
		{
            Point nZero = (f[1] <= -threshold) ? \
                           n1 : \
                           Point({ (f[1] * n0.x() - f[0] * n1.x())/ (f[1] - f[0]), (f[1] * n0.y() - f[0] * n1.y())/ (f[1] - f[0])});

            h = (nZero - n0).length();

            intRho = 0.5*h*(
                        cell.reconstructSolution(n0, 0)    - edgeBound.bc->applyBoundary(cell.reconstructSolution(n0))[0] +
                        cell.reconstructSolution(nZero, 0) - edgeBound.bc->applyBoundary(cell.reconstructSolution(nZero))[0]
                  );

            intE = 0.5*h*(
                        cell.reconstructSolution(n0, 4)    - edgeBound.bc->applyBoundary(cell.reconstructSolution(n0))[4] +
                        cell.reconstructSolution(nZero, 4) - edgeBound.bc->applyBoundary(cell.reconstructSolution(nZero))[4]
                  );

            intP = 0.5*h*(
                        problem.getPressure(cell.reconstructSolution(n0))    - problem.getPressure(edgeBound.bc->applyBoundary(cell.reconstructSolution(n0))) +
                        problem.getPressure(cell.reconstructSolution(nZero)) - problem.getPressure(edgeBound.bc->applyBoundary(cell.reconstructSolution(nZero)))
                  );


            if (debugFlux)
            {
//#pragma omp master
                cout << "\n----\nu0 < 0 bound \n -- \n";
                cout << "h = " << h << endl;
                cout << "nZero = " << nZero << endl;
                cout << "mySol zero rho = "   << cell.reconstructSolution(nZero, 0) << endl;
                cout << "neibSol zero rho = " << edgeBound.bc->applyBoundary(cell.reconstructSolution(nZero))[0] << endl;
                cout << "mySol zero e = "     << cell.reconstructSolution(nZero, 4) << endl;
                cout << "neibSol zero e = "   << edgeBound.bc->applyBoundary(cell.reconstructSolution(nZero))[4] << endl;
                cout << "mySol n0 rho = "     << cell.reconstructSolution(n0, 0) << endl;
                cout << "neibSol n0 rho = "   << edgeBound.bc->applyBoundary(cell.reconstructSolution(n0))[0] << endl;
                cout << "mySol n0 e = "       << cell.reconstructSolution(n0, 4) << endl;
                cout << "neibSol n0 e = "     << edgeBound.bc->applyBoundary(cell.reconstructSolution(n0))[4] << endl;

                cout << "int rho = "          << intRho << endl;
                cout << "int e = "            << intE << endl;
            }
        }
		else
		{
            Point nZero = Point({ (f[1] * n0.x() - f[0] * n1.x())/ (f[1] - f[0]), (f[1] * n0.y() - f[0] * n1.y())/ (f[1] - f[0])});

            h = (n1 - nZero).length();

            intRho = 0.5*h*(
                        cell.reconstructSolution(n1, 0)    - edgeBound.bc->applyBoundary(cell.reconstructSolution(n1))[0] +
                        cell.reconstructSolution(nZero, 0) - edgeBound.bc->applyBoundary(cell.reconstructSolution(nZero))[0]
                  );

            intE = 0.5*h*(
                        cell.reconstructSolution(n1, 4)    - edgeBound.bc->applyBoundary(cell.reconstructSolution(n1))[4] +
                        cell.reconstructSolution(nZero, 4) - edgeBound.bc->applyBoundary(cell.reconstructSolution(nZero))[4]
                  );

            intP = 0.5*h*(
                        problem.getPressure(cell.reconstructSolution(n1))    - problem.getPressure(edgeBound.bc->applyBoundary(cell.reconstructSolution(n1))) +
                        problem.getPressure(cell.reconstructSolution(nZero)) - problem.getPressure(edgeBound.bc->applyBoundary(cell.reconstructSolution(nZero)))
                  );

            if (debugFlux)
            {
//#pragma omp master
                cout << "\n----\nu0 > 0, u1 < 0 bound\n -- \n";

                cout << "h = " << h << endl;
                cout << "nZero = " << nZero << endl;
                cout << "mySol zero rho = "   << cell.reconstructSolution(nZero, 0) << endl;
                cout << "neibSol zero rho = " << edgeBound.bc->applyBoundary(cell.reconstructSolution(nZero))[0] << endl;
                cout << "mySol zero e = "     << cell.reconstructSolution(nZero, 4) << endl;
                cout << "neibSol zero e = "   << edgeBound.bc->applyBoundary(cell.reconstructSolution(nZero))[4] << endl;
                cout << "mySol n0 rho = "     << cell.reconstructSolution(n1, 0) << endl;
                cout << "neibSol n0 rho = "   << edgeBound.bc->applyBoundary(cell.reconstructSolution(n1))[0] << endl;
                cout << "mySol n0 e = "       << cell.reconstructSolution(n1, 4) << endl;
                cout << "neibSol n0 e = "     << edgeBound.bc->applyBoundary(cell.reconstructSolution(n1))[4] << endl;

                cout << "int rho = "          << intRho << endl;
                cout << "int e = "            << intE << endl;
            }
		}

        return
        {
            intRho,
            intE,
            intP,
            h
        };
    } // end catch
}

bool checkNegativeScalar(const Cell& cell)
{
    bool negScalar = false;
    
    for (int i = 0; i < cell.nEntities; ++i)
    {
        //const shared_ptr<Point> node = cell.nodes[i];
	numvector<double, 5> sol = cell.reconstructSolution(cell.nodes[i]);
        
        if (sol[0] < threshold || \
            sol[4] < threshold || \
            cell.problem.getPressure(sol) < threshold)
        {
            negScalar = true;
            break;
        }
    }
    
    return negScalar;
}


                                                

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


vector<int> IndicatorKXRCF::checkDiscontinuities() const
{
    vector<int> troubledCells;
    
    troubledCells.reserve(mesh.nCells);
    
    double indicatorRho = 0.0;
    double indicatorE = 0.0;
    double indicatorP = 0.0;
    bool negScalar = false;


//    cout << "\n========\nRun KXRCF\n--------\n";
#pragma omp parallel for \
shared(mesh, troubledCells, cout) \
firstprivate(indicatorRho, indicatorE, indicatorP, negScalar) \
default (none)
    for (int i = 0; i < mesh.nCells; ++i)
    //for (const shared_ptr<Cell> cell : mesh.cells)
    {
        
        const shared_ptr<Cell> cell = mesh.cells[i];
//        if (cell->number == 90 || cell->number == 148 || cell->number == 1026 || cell->number == 1468)
//            debugFlux = true;
//        else
//            debugFlux = false;

//        if (cell->number == 190 || cell->number == 206)
//            debugFlux = true;
//        else
//            debugFlux = false;

        // check negative rho/E values

        //for (const shared_ptr<Point> node : cell->nodes)
        
        negScalar = checkNegativeScalar(*cell);

        if (negScalar)
        {
            //cout << "true neg scalar on tthreasd " << omp_get_thread_num() << endl;
            negScalar = false;
            troubledCells.emplace_back(cell->number);
            
        }
        else
        {

        // get ~ radius of circumscribed circle in cell[i]
        double h = sqrt(0.5*cell->getArea());

		// get norm Q_j
        double normQrho = cell->getNormQ(0);
        double normQe = cell->getNormQ(4);
        double normQp = cell->getNormQp();

		// get momentum in cell nodes
        //numvector<double, 4> integral = {0.0, 0.0, 0.0, 0.0}; // rho|e|p|length
        numvector<double, 4> integral; // rho|e|p|length
        integral[0] = 0.0;
        integral[1] = 0.0;
        integral[2] = 0.0;
        integral[3] = 0.0;

        //for (const shared_ptr<Edge> edge : cell->edges0)
        for (int i = 0; i < cell->nEntities; ++i)
        {
            const shared_ptr<Edge> edge = cell->edges[i];
            integral += massFlux(*edge, *cell);
        }

        indicatorRho = fabs(integral[0]) / max((h*integral[3]*normQrho), threshold);
        indicatorE   = fabs(integral[1]) / max((h*integral[3]*normQe), threshold);
        indicatorP   = fabs(integral[2]) / max((h*integral[3]*normQp), threshold);

        if (debugFlux)
        {
//#pragma omp master
            cout << "integralRho = " << fabs(integral[0]) << ", integralE = " << fabs(integral[1]) << ", len = " << integral[2] << endl;
        }

        if (debugFlux)
        {
//#pragma omp master
            cout << "cell #" << cell->number <<": indicatorRho = " << indicatorRho << ", indicatorE = " << indicatorE << "\n===============\n";
        }

        if (indicatorRho > 1.0 || indicatorE > 1.0 || indicatorP > 1.0)
            troubledCells.emplace_back(cell->number);
	}
//        if (indicatorRho > 1.0)
//            cout << "cell #" << cell->number << ", ind = " << indicatorRho << endl;
    }

//    cout << "\ntroubled cells: " ;

//    for (int iCell : troubledCells)
//        cout << iCell << ' ';

//    cout << endl;


//    vector<int> isTroubled(mesh.nCells);

//    for (int iCell : troubledCells)
//        isTroubled[iCell] = 1;

//    ofstream writer;
//    writer.open("troubledCells");

//    writer << "SCALARS isTroubled int" << endl;
//    writer << "LOOKUP_TABLE default" << endl;

//    for(int i = 0; i < mesh.nCells; ++i)
//        writer << isTroubled[i] << endl;
    
//    writer.close();

    return troubledCells;
}
