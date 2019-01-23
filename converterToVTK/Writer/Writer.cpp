#include "Writer.h"
#include "Params.h"

using namespace std;

void Writer::exportMeshVTK(const std::string& fileName) const
{
    ofstream wStream (fileName);
    exportMeshVTK(wStream);
    wStream.close();
}


void Writer::exportFrameVTK(const std::string& fileName) const
{
    ofstream wStream (fileName);
    exportMeshVTK(wStream);
    exportFrameVTK(wStream);
    wStream.close();
}


void Writer::exportMeshVTK(ostream& wStream) const
{
    //wStream.open("Mesh.vtk");


    wStream << "# vtk DataFile Version 2.0" << endl;
    wStream << "RKDG 2D data"               << endl;
    wStream << "ASCII"                      << endl;

    wStream << "DATASET POLYDATA"           << endl;

    wStream << "POINTS " << mesh.nNodes << " float" << endl;

    // write coordinates of nodes
    for (int i = 0; i < mesh.nNodes; ++i)
        wStream << mesh.nodes[i]->x() << ' ' << mesh.nodes[i]->y() << ' ' << "0" << endl;

    //get size of polygon list

    int polySize = 0;

    for (int i = 0; i < mesh.nRealCells; ++i)
        polySize += mesh.cells[i]->nEntities;

    polySize += mesh.nRealCells;

    wStream << "POLYGONS " << mesh.nRealCells << ' ' << polySize << endl;


    // write cells using numbers of nodes
    for (int i = 0; i < mesh.nRealCells; ++i)
    {
        wStream << mesh.cells[i]->nEntities << ' ';

        for (const shared_ptr<Point> node : mesh.cells[i]->nodes)
            wStream << node->gNumber << ' ';

        wStream << endl;
    }

}


void Writer::exportFrameVTK(ostream& wStream) const
{
    // get cell data

    wStream << "CELL_DATA " << mesh.nRealCells << endl;

    wStream << "SCALARS rho double" << endl;
    wStream << "LOOKUP_TABLE default" << endl;

    //  for (const shared_ptr<Cell> cell : mesh.cells)
    for (int iCell = 0; iCell < mesh.nRealCells; ++iCell)
        wStream << solution.reconstruct(iCell, mesh.cells[iCell]->getCellCenter(), RHO) << endl;

    wStream << "SCALARS e double" << endl;
    wStream << "LOOKUP_TABLE default" << endl;

    for (int iCell = 0; iCell < mesh.nRealCells; ++iCell)
        wStream << solution.reconstruct(iCell, mesh.cells[iCell]->getCellCenter(), E) << endl;

    wStream << "SCALARS p double" << endl;
    wStream << "LOOKUP_TABLE default" << endl;

    for (int iCell = 0; iCell < mesh.nRealCells; ++iCell)
        wStream << physics.getPressure(solution.reconstruct(iCell, mesh.cells[iCell]->getCellCenter())) << endl;


    wStream << "VECTORS U double" << endl;

    for (int iCell = 0; iCell < mesh.nRealCells; ++iCell)
    {
        double rho = solution.reconstruct(iCell, mesh.cells[iCell]->getCellCenter(), RHO);
        wStream << solution.reconstruct(iCell, mesh.cells[iCell]->getCellCenter(), RHOU) / rho << endl;
        wStream << solution.reconstruct(iCell, mesh.cells[iCell]->getCellCenter(), RHOV) / rho << endl;
        wStream << 0.0 << endl;
    }

    //indicator.writeTroubledCellsVTK();

    // get point data

    int nEntitiesTotal = 0;

    for (int i = 0; i < mesh.nRealCells; ++i)
        nEntitiesTotal += mesh.cells[i]->nEntities;

    wStream << "POINT_DATA " << nEntitiesTotal << endl;
    
    /*
    wStream << "SCALARS rho double" << endl;
    wStream << "LOOKUP_TABLE default" << endl;

    for (const shared_ptr<Cell> cell : mesh.cells)
    for (int j = 0; j < cell->nEntities; ++j)
    wStream << cell->reconstruct(cell->nodes[j], 0) << endl;

    wStream << "SCALARS e double" << endl;
    wStream << "LOOKUP_TABLE default" << endl;

    for (const shared_ptr<Cell> cell : mesh.cells)
    for (int j = 0; j < cell->nEntities; ++j)
    wStream << cell->reconstruct(cell->nodes[j], 4) << endl;

    wStream << "SCALARS p double" << endl;
    wStream << "LOOKUP_TABLE default" << endl;

    for (const shared_ptr<Cell> cell : mesh.cells)
    for (int j = 0; j < cell->nEntities; ++j)
    wStream << problem.getPressure(cell->reconstruct(cell->nodes[j])) << endl;


    wStream << "VECTORS U double" << endl;

    for (const shared_ptr<Cell> cell : mesh.cells)
    for (int j = 0; j < cell->nEntities; ++j)
    {
    double rho = cell->reconstruct(cell->nodes[j], 0);
    wStream << cell->reconstruct(cell->nodes[j], 1) / rho << endl;
    wStream << cell->reconstruct(cell->nodes[j], 2) / rho << endl;
    wStream << 0.0 << endl;
    }
    */
}
