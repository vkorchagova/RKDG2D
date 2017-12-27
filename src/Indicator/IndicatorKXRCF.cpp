#include "IndicatorKXRCF.h"

using namespace std;

Point massFlux(const Edge& edge, const Cell& cell, const Point& nrm)
{
    double mom[2];
    Point tau;
        
	try
	{
        // for internal edges
        const EdgeInternal& edgeInt = dynamic_cast<const EdgeInternal&>(edge);

        // neighbour cell could not be equal to considered cell
        Cell& neib = (&cell == edge.neibCells[0].get()) ? *edge.neibCells[1] : *edge.neibCells[0];

		if (fabs(nrm.x()) < 1e-10)
        { // if edge is hor
			mom[0] = cell.reconstructSolution(edge.nodes[0], 2) * nrm.y();
			mom[1] = cell.reconstructSolution(edge.nodes[1], 2) * nrm.y();

			tau[0] = 1.0; tau[1] = 0.0;
		}
		else
        { // now -- if edge is ver
			mom[0] = cell.reconstructSolution(edge.nodes[0], 1) * nrm.x();
			mom[1] = cell.reconstructSolution(edge.nodes[1], 1) * nrm.x();

			tau[0] = 0.0; tau[1] = 1.0;
		}


        // 4 variants of flux integral (only through the part of edge with inward flux)
        // we have linear functions -> use trapezia formulae for integration

        if ((mom[0] >= 0.0) && (mom[1] >= 0.0))
			return Point({ 0.0, 0.0 });

		else if (mom[0] < 0.0)
		{
            double h = (mom[1] <= 0.0) ? edge.getLength() : -(mom[0] * edge.getLength() / (mom[1] - mom[0]));

            double mySolH = cell.reconstructSolution(*edge.nodes[0] + tau * h, 0);
            double neibSolH = neib.reconstructSolution(*edge.nodes[0] + tau * h, 0);

			return Point(
			{ 0.5*h*(cell.reconstructSolution(edge.nodes[0], 0) - neib.reconstructSolution(edge.nodes[0], 0) + mySolH - neibSolH), h }
			);
		}
		else
		{
			double h = -(mom[0] * edge.getLength() / (mom[1] - mom[0]));
			double mySolH = cell.reconstructSolution(*edge.nodes[0] + tau*h, 0);
			double neibSolH = neib.reconstructSolution(*edge.nodes[0] + tau*h, 0);

			return Point(
			{ 0.5*(edge.getLength() - h)*(cell.reconstructSolution(edge.nodes[1], 0) - neib.reconstructSolution(edge.nodes[1], 0) + mySolH - neibSolH), edge.getLength() - h }
			);
		}

		cout << "STRANGE!!!" << endl;
		return Point({ 0.0, 0.0 });

	}
	catch (...)
	{
        // for edge on boundary
        const EdgeBoundary& edgeBound = dynamic_cast<const EdgeBoundary&>(edge);

		//Cell& neib = (&cell == edge.neibCells[0].get()) ? *edge.neibCells[1] : *edge.neibCells[0];
		if (fabs(nrm.x()) < 1e-10)
		{
			mom[0] = cell.reconstructSolution(edge.nodes[0], 2) * nrm.y();
			mom[1] = cell.reconstructSolution(edge.nodes[1], 2) * nrm.y();

			tau[0] = 1.0; tau[1] = 0.0;
		}
		else
		{
			mom[0] = cell.reconstructSolution(edge.nodes[0], 1) * nrm.x();
			mom[1] = cell.reconstructSolution(edge.nodes[1], 1) * nrm.x();

			tau[0] = 0.0; tau[1] = 1.0;
		}

		if ((mom[0] >= 0.0) && (mom[1] >= 0.0))
			return Point({ 0.0, 0.0 });

		else if (mom[0] < 0.0)
		{
            double h = (mom[1] <= 0) ? edge.getLength() : -(mom[0] * edge.getLength() / (mom[1] - mom[0]));
			double mySolH = cell.reconstructSolution(*edge.nodes[0] + tau*h, 0);
			double neibSolH = edgeBound.applyBoundary(cell.reconstructSolution(*edge.nodes[0] + tau*h))[0];// neib.reconstructSolution( *edge.nodes[0] + tau*h, 0);

			return Point({ 0.5*h*(cell.reconstructSolution(edge.nodes[0], 0) - edgeBound.applyBoundary(cell.reconstructSolution(*edge.nodes[0]))[0] + mySolH - neibSolH), h });
		}
		else
		{
			double h = -(mom[0] * edge.getLength() / (mom[1] - mom[0]));
			double mySolH = cell.reconstructSolution(*edge.nodes[0] + tau*h, 0);
			double neibSolH = edgeBound.applyBoundary(cell.reconstructSolution(*edge.nodes[0] + tau*h))[0];// neib.reconstructSolution( *edge.nodes[0] + tau*h, 0);

			return Point(
			{ 0.5*(edge.getLength() - h)*(cell.reconstructSolution(edge.nodes[1], 0) - edgeBound.applyBoundary(cell.reconstructSolution(*edge.nodes[1]))[0] + mySolH - neibSolH), edge.getLength() - h }
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


vector<double> IndicatorKXRCF::checkDiscontinuities() const
{
    vector<double> indicator(mesh.nCells);
    
	for (int i = 0; i < mesh.nCells; ++i)
	{
		// get ~ radius of circumscribed circle in cell[i]
		double h = sqrt(0.5 * mesh.cells[i]->getArea());

		// get norm Q_j
		double normQ = mesh.cells[i]->getNormQ();

		// get momentum in cell nodes

        //vector<shared_ptr<Point>> nodes = mesh.cells[i]->getCellCoordinates(); // store nodes?..

        //vector<Point> momentum(nodes.size());

		const Cell& mci = *mesh.cells[i];


        //for (size_t j = 0; j < nodes.size(); ++j)
        //{
        //	momentum[j] = Point({ mci.reconstructSolution(nodes[j], 1), mci.reconstructSolution(nodes[j], 2) });
        //}

		Point integral(0.0);

		// for bottom edge
        integral += massFlux(*mci.edges[0], mci, Point({ 0.0, -1.0 }));

		// for top edge
		integral += massFlux(*mci.edges[2], mci, Point({ 0.0, 1.0 }));

		// for right edge
        integral += massFlux(*mci.edges[1], mci, Point({ 1.0, 0.0 }));

		// for left edge
        integral += massFlux(*mci.edges[3], mci, Point({ -1.0, 0.0 }));


		indicator[i] = fabs(integral.x()) / max((h*integral.y()*normQ), 1e-10);

	}
    
    return indicator;
}
