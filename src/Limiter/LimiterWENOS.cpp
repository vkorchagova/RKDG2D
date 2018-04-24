#include "LimiterWENOS.h"

using namespace std;

void LimiterWENOS::limit(vector<numvector<double, 5 * nShapes>>& alpha)
{
    return limitTog(alpha);
}

void LimiterWENOS::limitTog(vector<numvector<double, 5 * nShapes>>& alpha)
{
    problem.setAlpha(alpha);

    vector<int> troubledCells = indicator.checkDiscontinuities();

    vector<double> gamma;
    double g = 0.001;

    // smoothness indicators

    vector<numvector<double, 5>> beta;

    // nonlinear weights

    vector<numvector<double, 5>> w;
    vector<numvector<double, 5>> wTilde;
    numvector<double, 5> wSum;

    // mean values

    vector<numvector<double, 5>> uMean;

    // p polynoms

    vector<numvector<double, 5 * nShapes>> p;


    // limit solution in troubled cells

    for (int iCell : troubledCells)
    {


        // find neighbours

        shared_ptr<Cell> cell = indicator.mesh.cells[iCell];

        vector<shared_ptr<Cell>> neibCellsX = cell->findNeighbourCellsX();
        vector<shared_ptr<Cell>> neibCellsY = cell->findNeighbourCellsY();

        //cout << "coeffs right neighbour = " << alpha[neibCellsX[1]->number] << endl;

        vector<shared_ptr<Cell>> cells = { cell };


        cells.insert(cells.end(), neibCellsX.begin(), neibCellsX.end());
        cells.insert(cells.end(), neibCellsY.begin(), neibCellsY.end());

        int nCells = cells.size();

        // get mean values

        uMean.resize(nCells);

        for (size_t k = 0; k < nCells; ++k)
            uMean[k] = cells[k]->reconstructSolution(cells[0]->getCellCenter());

        // get coeffs for polynoms p

        p.resize(nCells);

        for (size_t k = 0; k < nCells; ++k)
            p[k] = alpha[cells[k]->number] ;

        for (size_t k = 0; k < nCells; ++k)
            for (int i = 0; i < 5; ++i)
                for (int j = 1; j < nShapes; ++j)
                    p[k][i*nShapes + j] *= cells[k]->offsetPhi[j];

        for (size_t k = 0; k < nCells; ++k)
            for (int i = 0; i < 5; ++i)
                p[k][i*nShapes] = - uMean[k][i] + uMean[0][i] + cells[k]->reconstructSolution(cells[k]->getCellCenter(),i);

        // get linear weights

        gamma.resize(nCells);

        gamma[0] = 1.0 - (nCells - 1) * g;

        for (size_t k = 1; k < nCells; ++k)
            gamma[k] = g;

        // get smoothness indicators and nonlinear weights

        beta.resize(nCells);
        wTilde.resize(nCells);
        w.resize(nCells);

//        for (size_t k = 0; k < nCells; ++k)
//            for (int j = 0; j < 5; ++j)
//            {
//                beta[k][j] =  (cells[0]->h().x() * cells[0]-> h().y()) *(sqr(p[k][j*nShapes + 1]) + sqr(p[k][j*nShapes + 2]));
//                wTilde[k][j] = gamma[k] * (1.0 / sqr(beta[k][j] + 1e-6));
//            }



//        wSum = {0.0, 0.0, 0.0, 0.0, 0.0};

//        for (int j = 0; j < 5; ++j)
//            for (size_t k = 0; k < nCells; ++k)
//                wSum[j] += wTilde[k][j];

//        //cout << wSum << endl;

//        for (size_t k = 0; k < nCells; ++k)
//            for (int j = 0; j < 5; ++j)
//                w[k][j] = wTilde[k][j] / wSum[j];

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


        // project limited solution onto cell basis

        function<numvector<double, 5>(const Point& r)> foo = [=](const Point& r) \
        {
            numvector<double, 5> sum (0.0);

            for (int i = 0; i < 5; ++i)
                for (size_t k = 0; k < nCells; ++k)
                    sum[i] += w[k][i] * (p[k][i*nShapes] + \
                                         p[k][i*nShapes + 1] * (r.x() - cells[k]->getCellCenter().x()) + \
                                         p[k][i*nShapes + 2] * (r.y() - cells[k]->getCellCenter().y()));


            return sum;
        };

        alpha[iCell] = cells[0]->projection(foo);

        problem.setAlpha(alpha);

    }

}

void LimiterWENOS::limitSep(vector<numvector<double, 5 * nShapes>>& alpha)
{
    problem.setAlpha(alpha);
    vector<int> troubledCells = indicator.checkDiscontinuities();

    // limit in one-dimensional dir

    numvector<double, 5 * nShapes> alphax;
    numvector<double, 5 * nShapes> alphay;

    // linear weights

    vector<double> gammax;
    vector<double> gammay;
    double g = 0.001;

    // smoothness indicators

    vector<numvector<double, 5>> betax;
    vector<numvector<double, 5>> betay;

    // nonlinear weights

    vector<numvector<double, 5>> wx;
    vector<numvector<double, 5>> wy;
    vector<numvector<double, 5>> wTildex;
    vector<numvector<double, 5>> wTildey;
    numvector<double, 5> wSumx;
    numvector<double, 5> wSumy;

    // mean values

    vector<numvector<double, 5>> uMeanX;
    vector<numvector<double, 5>> uMeanY;

    // p polynoms

    vector<numvector<double, 5 * nShapes>> px;
    vector<numvector<double, 5 * nShapes>> py;


    // limit solution in troubled cells

    for (int iCell : troubledCells)
    {
        shared_ptr<Cell> cell = indicator.mesh.cells[iCell];

        vector<shared_ptr<Cell>> neibCellsX = cell->findNeighbourCellsX();
        vector<shared_ptr<Cell>> neibCellsY = cell->findNeighbourCellsY();

        vector<shared_ptr<Cell>> cellsHor = { cell };
        vector<shared_ptr<Cell>> cellsVer = { cell };

        cellsHor.insert(cellsHor.end(), neibCellsX.begin(), neibCellsX.end());
        cellsVer.insert(cellsVer.end(), neibCellsY.begin(), neibCellsY.end());

        size_t nCellsHor = cellsHor.size();
        size_t nCellsVer = cellsVer.size();

        // get mean values

        uMeanX.resize(nCellsHor);
        uMeanY.resize(nCellsVer);

        for (size_t i = 0; i < nCellsHor; ++i)
            uMeanX[i] = cellsHor[i]->reconstructSolution(cellsHor[0]->getCellCenter());

        for (size_t i = 0; i < cellsVer.size(); ++i)
            uMeanY[i] = cellsVer[i]->reconstructSolution(cellsVer[0]->getCellCenter());

        // get coeffs for polynoms p

        px.resize(nCellsHor);
        py.resize(nCellsVer);

        for (size_t i = 0; i < cellsHor.size(); ++i)
            px[i] = alpha[cellsHor[i]->number];

        for (size_t i = 0; i < cellsVer.size(); ++i)
            py[i] = alpha[cellsVer[i]->number];

        for (size_t k = 0; k < nCellsHor; ++k)
            for (int i = 0; i < 5; ++i)
                for (int j = 0; j < nShapes; ++j)
                    px[k][i*nShapes + j] *= cellsHor[k]->offsetPhi[j];

        for (size_t k = 0; k < nCellsVer; ++k)
            for (int i = 0; i < 5; ++i)
                for (int j = 0; j < nShapes; ++j)
                    py[k][i*nShapes + j] *= cellsVer[k]->offsetPhi[j];

        for (size_t k = 0; k < nCellsHor; ++k)
            for (int i = 0; i < 5; ++i)
                px[k][i*nShapes] = - uMeanX[k][i] + uMeanX[0][i] + cellsHor[k]->reconstructSolution(cellsHor[k]->getCellCenter(),i);

        for (size_t k = 0; k < nCellsVer; ++k)
            for (int i = 0; i < 5; ++i)
                py[k][i*nShapes] = - uMeanY[k][i] + uMeanY[0][i] + cellsVer[k]->reconstructSolution(cellsVer[k]->getCellCenter(),i);

        // get linear weights
        gammax.resize(nCellsHor);
        gammay.resize(nCellsVer);

        gammax[0] = 1.0 - (nCellsHor - 1) * g;
        gammay[0] = 1.0 - (nCellsVer - 1) * g;

        for (size_t k = 1; k < nCellsHor; ++k)
            gammax[k] = g;

        for (size_t k = 1; k < nCellsVer; ++k)
            gammay[k] = g;

        // get smoothness indicators and nonlinear weights

        betax.resize(nCellsHor);
        wTildex.resize(nCellsHor);
        wx.resize(nCellsHor);

        betay.resize(nCellsVer);
        wTildey.resize(nCellsVer);
        wy.resize(nCellsVer);

//        for (size_t k = 0; k < nCellsHor; ++k)
//        {
//            for (int j = 0; j < 5; ++j)
//            {
//                betax[k][j] = cellsHor[k]->h().x() * cellsHor[0]->h().x() * sqr(px[k][j*nShapes + 1]);
//                wTildex[k][j] = gammax[k] * (1.0 / sqr(betax[k][j] + 1e-6));
//            }
//        }

//        for (size_t k = 0; k < nCellsVer; ++k)
//        {
//            for (int j = 0; j < 5; ++j)
//            {
//                betay[k][j] = cellsVer[k]->h().y() * cellsVer[0]->h().y() * sqr( py[k][j*nShapes + 2]);
//                wTildey[k][j] = gammay[k] * (1.0 / sqr(betay[k][j] + 1e-6));
//            }
//        }

        wSumx = { 0.0, 0.0, 0.0, 0.0, 0.0 };
        wSumy = { 0.0, 0.0, 0.0, 0.0, 0.0 };

        for (size_t k = 0; k < nCellsHor; ++k)
            wSumx += wTildex[k];

        for (size_t k = 0; k < nCellsVer; ++k)
            wSumy += wTildey[k];

        for (size_t k = 0; k < nCellsHor; ++k)
            for (int j = 0; j < 5; ++j)
                wx[k][j] = wTildex[k][j] / wSumx[j];

        for (size_t k = 0; k < nCellsVer; ++k)
            for (int j = 0; j < 5; ++j)
                wy[k][j] = wTildey[k][j] / wSumy[j];

        // project limited solution onto cell basis

        function<numvector<double, 5>(const Point& r)> foox = [=](const Point& r) \
        {
            numvector<double, 5> sum (0.0);

            for (int i = 0; i < 5; ++i)
                for (size_t k = 0; k < nCellsHor; ++k)
                    sum[i] += wx[k][i] * (px[k][i*nShapes] + \
                                         px[k][i*nShapes + 1] * (r.x() - cellsHor[k]->getCellCenter().x()) + \
                                         px[k][i*nShapes + 2] * (r.y() - cellsHor[k]->getCellCenter().y()));


            return sum;
        };

        function<numvector<double, 5>(const Point& r)> fooy = [=](const Point& r) \
        {
            numvector<double, 5> sum (0.0);

            for (int i = 0; i < 5; ++i)
                for (size_t k = 0; k < nCellsVer; ++k)
                    sum[i] += wy[k][i] * (py[k][i*nShapes] + \
                                         py[k][i*nShapes + 1] * (r.x() - cellsVer[k]->getCellCenter().x()) + \
                                         py[k][i*nShapes + 2] * (r.y() - cellsVer[k]->getCellCenter().y()));


            return sum;
        };

        alphax = cellsHor[0]->projection(foox);
        alphay = cellsVer[0]->projection(fooy);

        alpha[iCell] = (alphax + alphay) * 0.5;
        //alpha[iCell] = alphax;

        problem.setAlpha(alpha);
    }

    return;
}

void LimiterWENOS::limitX(vector<numvector<double, 5 * nShapes>>& alpha)
{
    problem.setAlpha(alpha);
    vector<int> troubledCells = indicator.checkDiscontinuities();

    // linear weights

    numvector<double, 3> gamma = { 0.001, 0.998, 0.001 };

    // smoothness indicators

    numvector<numvector<double, 5>, 3> beta;

    // nonlinear weights

    numvector<numvector<double, 5>, 3> w;
    numvector<numvector<double, 5>, 3> wTilde;
    numvector<double, 5> wSum;

    // mean values

    numvector<numvector<double, 5>, 3> uMean;

    // p polynoms

    numvector<numvector<double, 5 * nShapes>, 3> p;


    // limit solution in troubled cells

    for (int icell : troubledCells)
    {
        problem.setAlpha(alpha);
        // limit solution in X direction

        // find neighbours in x direction (we know the mesh is rectilinear)

        numvector<shared_ptr<Cell>, 3> cellsHor = { indicator.mesh.cells[icell-1], \
                                                    indicator.mesh.cells[icell], \
                                                    indicator.mesh.cells[icell+1] };

        // get mean values

        for (int k = 0; k < 3; ++k)
            uMean[k] = cellsHor[k]->reconstructSolution(cellsHor[1]->getCellCenter());

        // get coeffs for polynoms p

        p[0] = alpha[icell-1];
        p[1] = alpha[icell];
        p[2] = alpha[icell+1];

        for (int k = 0; k < 3; ++k)
            for (int i = 0; i < 5; ++i)
                for (int j = 0; j < nShapes; ++j)
                    p[k][i*nShapes + j] *= cellsHor[k]->offsetPhi[j];

        for (int i = 0; i < 5; ++i)
            p[0][i*nShapes] = - uMean[0][i] + uMean[1][i] + cellsHor[0]->reconstructSolution(cellsHor[0]->getCellCenter(),i);

        for (int i = 0; i < 5; ++i)
            p[2][i*nShapes] = - uMean[2][i] + uMean[1][i] + cellsHor[2]->reconstructSolution(cellsHor[2]->getCellCenter(),i);

        // get smoothness indicators and nonlinear weights

//        for (int k = 0; k < 3; ++k)
//        {
//            for (int j = 0; j < 5; ++j)
//            {
//                beta[k][j] = sqr( min(cellsHor[k]->h().x(),cellsHor[k]->h().y()) * p[k][j*nShapes + 1] );
//                wTilde[k][j] = gamma[k] * (1.0 / sqr(beta[k][j] + 1e-6));
//            }
//        }

//        wSum = wTilde[0] + wTilde[1] + wTilde[2];

//        for (int k = 0; k < 3; ++k)
//        {
//            for (int j = 0; j < 5; ++j)
//                w[k][j] = wTilde[k][j] / wSum[j];
//        }

        // project limited solution onto cell basis

        function<numvector<double, 5>(const Point& r)> foo = [=](const Point& r) \
        {
            numvector<double, 5> sum (0.0);

            for (int i = 0; i < 5; ++i)
                for (int k = 0; k < 3; ++k)
                    sum[i] += w[k][i] * (p[k][i*nShapes] + \
                                         p[k][i*nShapes + 1] * (r.x() - cellsHor[k]->getCellCenter().x()) + \
                                         p[k][i*nShapes + 2] * (r.y() - cellsHor[k]->getCellCenter().y()));


            return sum;
        };

        alpha[icell] = cellsHor[1]->projection(foo);
    }

    return;
}
