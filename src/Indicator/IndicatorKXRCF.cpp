#include "IndicatorKXRCF.h"

using namespace std;

bool debugFlux;

Point massFlux(const Edge& edge, const Cell& cell)
{
    double f[2];

    Point n = Point({0.0, 0.0});

    Point n0 = *edge.nodes[0];
    Point n1 = *edge.nodes[1];

    double threshold = 1e-10;
        
	try
	{
        // for internal edges
        const EdgeInternal& edgeInt = dynamic_cast<const EdgeInternal&>(edge);

        // neighbour cell could not be equal to considered cell
        Cell& neib = (&cell == edge.neibCells[0].get()) ? *edge.neibCells[1] : *edge.neibCells[0];

        // correct normal position: outside related to considered cell
        n = ((n0 - cell.getCellCenter())*n > 0) ? edge.n : Point(-edge.n);

        //get normal velocity component in nodes of edge
        f[0] = rotate(cell.reconstructSolution(n0), n)[1];
        f[1] = rotate(cell.reconstructSolution(n1), n)[1];

        // 4 variants of flux integral (only through the part of edge with inward flux)
        // we have linear functions -> use trapezia formulae for integration


        // no entrance
        if ((f[0] >= -threshold) && (f[1] >= -threshold))
			return Point({ 0.0, 0.0 });

        //entrance at n[0]
        else if (f[0] < -threshold)
		{
            Point nZero = (f[1] <= -threshold) ? \
                           *edge.nodes[1] : \
                           Point({ (f[1] * n0.x() - f[0] * n1.x())/ (f[1] - f[0]), (f[1] * n0.y() - f[0] * n1.y())/ (f[1] - f[0])});

            double h = (nZero - n0).length();

            double mySolH   = cell.reconstructSolution(nZero, 0);
            double neibSolH = neib.reconstructSolution(nZero, 0);

            if (debugFlux)
            {
                cout << "h = " << h << endl;
                cout << "mySolH = " << mySolH << endl;
                cout << "neibSolH = " << neibSolH << endl;
            }
			return Point(
            { 0.5*h*(cell.reconstructSolution(n0, 0) - neib.reconstructSolution(n0, 0) + mySolH - neibSolH), h }
			);
		}
        else // no entrance at n0, entrance at n1
		{
            Point nZero = Point({ (f[1] * n0.x() - f[0] * n1.x())/ (f[1] - f[0]), (f[1] * n0.y() - f[0] * n1.y())/ (f[1] - f[0])});

            double h = (n1 - nZero).length();

            double mySolH = cell.reconstructSolution(nZero, 0);
            double neibSolH = neib.reconstructSolution(nZero, 0);

			return Point(
            { 0.5*h*(cell.reconstructSolution(n1, 0) - neib.reconstructSolution(n1, 0) + mySolH - neibSolH), h }
			);
		}

		cout << "STRANGE!!!" << endl;
		return Point({ 0.0, 0.0 });

	}
	catch (...)
	{
        // for edge on boundary
        const EdgeBoundary& edgeBound = dynamic_cast<const EdgeBoundary&>(edge);

        n = edge.n;

        //get normal velocity component in nodes of edge
        f[0] = rotate(cell.reconstructSolution(edge.nodes[0], 1), n)[1];
        f[1] = rotate(cell.reconstructSolution(edge.nodes[1], 1), n)[1];


        if ((f[0] >= -threshold) && (f[1] >= -threshold))
			return Point({ 0.0, 0.0 });

        else if (f[0] < -threshold)
		{
            Point nZero = (f[1] <= -threshold) ? \
                           *edge.nodes[1] : \
                           Point({ (f[1] * n0.x() - f[0] * n1.x())/ (f[1] - f[0]), (f[1] * n0.y() - f[0] * n1.y())/ (f[1] - f[0])});

            double h = (nZero - n0).length();

            double mySolH = cell.reconstructSolution(nZero, 0);
            double neibSolH = edgeBound.bc->applyBoundary(cell.reconstructSolution(nZero))[0];// neib.reconstructSolution( *edge.nodes[0] + tau*h, 0);

            return Point({ 0.5*h*(cell.reconstructSolution(n0, 0) - edgeBound.bc->applyBoundary(cell.reconstructSolution(n0))[0] + mySolH - neibSolH), h });
		}
		else
		{
            Point nZero = Point({ (f[1] * n0.x() - f[0] * n1.x())/ (f[1] - f[0]), (f[1] * n0.y() - f[0] * n1.y())/ (f[1] - f[0])});

            double h = (n1 - nZero).length();

            double mySolH = cell.reconstructSolution(nZero, 0);
            double neibSolH = edgeBound.bc->applyBoundary(cell.reconstructSolution(nZero))[0];// neib.reconstructSolution( *edge.nodes[0] + tau*h, 0);

			return Point(
            { 0.5*h*(cell.reconstructSolution(n1, 0) - edgeBound.bc->applyBoundary(cell.reconstructSolution(n1))[0] + mySolH - neibSolH), h }
			);
		}

		cout << "STRANGE!!!" << endl;
		return Point({ 0.0, 0.0 });
		

		//cout << "CATCH!!!" << endl;
		//return Point({ 100.0, 100.0 });
	}

    cout << "REALLY STRANGE!!!" << endl;
    return Point({ 0.0, 0.0 });
    
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


vector<int> IndicatorKXRCF::checkDiscontinuities() const
{
    vector<int> troubledCells;

    double indicator;
    
    for (const shared_ptr<Cell> cell : mesh.cells)
	{
		// get ~ radius of circumscribed circle in cell[i]
        double h = sqrt(0.5 * cell->getArea());

		// get norm Q_j
        double normQ = cell->getNormQ();

		// get momentum in cell nodes

        double integralRho = 0.0;
        double integralE = 0.0;
        double lenEntrance = 0.0;

        Point integral = Point({0.0, 0.0});

        for (const shared_ptr<Edge> edge : cell->edges)
        {
            integral += massFlux(*edge, *cell);
        }


        if (cell->number == 6 || cell->number == 7 || cell->number == 8)
            debugFlux = true;
        else
            debugFlux = false;


        indicator = fabs(integral.x()) / max((h*integral.y()*normQ), 1e-10);

        if (debugFlux)
        {
            cout << "integral = " << fabs(integral.x()) << ", len = " << integral.y() << endl;

        }

        cout << "cell #" << cell->number <<": indicator = " << indicator << "\n===============\n";

        if (indicator > 1.0)
            troubledCells.push_back(cell->number);
	}

    cout << "\ntroubled cells: " ;

    for (int iCell : troubledCells)
        cout << iCell << ' ';

    cout << endl;

    return troubledCells;
}
