#include "Writer.h"
#include "Params.h"

using namespace std;

void Writer::exportMeshVTK(const std::string& fileName) const
{
    ofstream wStream (fileName);
    exportMeshVTK(wStream);
    wStream.close();
}

void Writer::exportMeshVTK_polyvertices(const std::string& fileName) const
{
    ofstream wStream (fileName);
    exportMeshVTK_polyvertices(wStream);
    wStream.close();
}


void Writer::exportFrameVTK(const std::string& fileName) const
{
    cout << "Export frame " << fileName << "..." << endl;
    
    ofstream wStream (fileName);
    exportMeshVTK(wStream);
    exportFrameVTK(wStream);
    wStream.close();
    
    cout << endl;
}

void Writer::exportFrameVTK_polyvertices(const std::string& fileName) const
{
    cout << "Export frame " << fileName << "..." << endl;
    
    ofstream wStream (fileName);
    exportMeshVTK_polyvertices(wStream);
    exportFrameVTK_polyvertices(wStream);
    wStream.close();
    
    cout << endl;
}

void Writer::exportNativeCoeffs(const std::string& fileName) const
{
    ofstream wStream (fileName);
    outputFullNativeCoeffs(wStream);
    wStream.close();
}


void Writer::exportMeshVTK(ostream& wStream) const
{
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

    wStream.precision(16);

    wStream << "CELL_DATA " << mesh.nRealCells << endl;

    wStream << "SCALARS rho double" << endl;
    wStream << "LOOKUP_TABLE default" << endl;

    //  for (int i = 0; i < mesh.nRealCells; ++i)
    for (int iCell = 0; iCell < mesh.nRealCells; ++iCell)
        wStream << solution.solToExport[iCell][0] << endl;

    wStream << "SCALARS e double" << endl;
    wStream << "LOOKUP_TABLE default" << endl;

    for (int iCell = 0; iCell < mesh.nRealCells; ++iCell)
        wStream << solution.solToExport[iCell][4] << endl;

    wStream << "SCALARS p double" << endl;
    wStream << "LOOKUP_TABLE default" << endl;

    for (int iCell = 0; iCell < mesh.nRealCells; ++iCell)
        wStream << solution.solToExport[iCell][5] << endl;


    wStream << "VECTORS U double" << endl;

    for (int iCell = 0; iCell < mesh.nRealCells; ++iCell)
    {
        wStream << solution.solToExport[iCell][1] << endl;
        wStream << solution.solToExport[iCell][2] << endl;
        wStream << 0.0 << endl;
    }

    //indicator.writeTroubledCellsVTK();

    // get point data

    int nEntitiesTotal = 0;

    for (int i = 0; i < mesh.nRealCells; ++i)
        nEntitiesTotal += mesh.cells[i]->nEntities;

    wStream << "POINT_DATA " << nEntitiesTotal << endl;
    
}


void Writer::outputNativeCoeffs(ostream& wStream) const
{
    wStream.precision(16);

    for (size_t k = 0; k < solution.SOL.size(); ++k)
    {
        for (int i = 0; i < dimS; ++i)
            wStream << solution.SOL[k][i] << ' ';

        wStream << endl;
    }
}

void Writer::outputFullNativeCoeffs(ostream& wStream) const
{
    wStream.precision(16);
    
    for (size_t k = 0; k < solution.fullSOL.size(); ++k)
    {
        for (int i = 0; i < dimS; ++i)
            wStream << solution.fullSOL[k][i] << ' ';

        wStream << endl;
    }
}

void Writer::outputFullNativeCoeffs(std::ostream& wStream, std::vector<numvector<double, dimS>> coeffs) const
{
    wStream.precision(16);
    
    for (size_t k = 0; k < coeffs.size(); ++k)
    {
        for (int i = 0; i < dimS; ++i)
            wStream << coeffs[k][i] << ' ';

        wStream << endl;
    }
}


void Writer::exportMeshVTK_polyvertices(ostream& wStream) const
{
    //wStream.open("Mesh.vtk");

    wStream << "# vtk DataFile Version 2.0" << endl;
    wStream << "RKDG 2D data" << endl;
    wStream << "ASCII" << endl;

    wStream << "DATASET POLYDATA" << endl;

    int polySize = 0;

    for (int i = 0; i < mesh.nRealCells; ++i)
        polySize += mesh.cells[i]->nEntities;

    wStream << "POINTS " << polySize << " float" << endl;

    // write coordinates of nodes
    for (int i = 0; i < mesh.nRealCells; ++i)
        for (int j = 0; j < mesh.cells[i]->nEntities; ++j)
            wStream << mesh.cells[i]->nodes[j]->x() << ' ' << mesh.cells[i]->nodes[j]->y() << ' ' << "0" << endl;

    //get size of polygon list

    polySize += mesh.nRealCells;

    wStream << "POLYGONS " << mesh.nRealCells << ' ' << polySize << endl;


    // write cells using numbers of nodes
    int numPolyVertex = 0;
    for (int i = 0; i < mesh.nRealCells; ++i)
    {
        wStream << mesh.cells[i]->nEntities << ' ';

        for (const shared_ptr<Point> node : mesh.cells[i]->nodes)
        {
            wStream << numPolyVertex << ' ';
            numPolyVertex++;
        }

        wStream << endl;
    }
}


void Writer::exportFrameVTK_polyvertices(ostream& wStream) const
{
    int polySize = 0;

    for (int i = 0; i < mesh.nRealCells; ++i)
        polySize += mesh.cells[i]->nEntities;

    // get cell data

    wStream << "CELL_DATA " << mesh.nRealCells << endl;

    wStream << "SCALARS rho double" << endl;
    wStream << "LOOKUP_TABLE default" << endl;

    for (int i = 0; i < mesh.nRealCells; ++i)
        wStream << solution.reconstruct(mesh.globalCellNumber[i], mesh.cells[i]->getCellCenter(), RHO) << endl;

    wStream << "SCALARS e double" << endl;
    wStream << "LOOKUP_TABLE default" << endl;

    for (int i = 0; i < mesh.nRealCells; ++i)
        wStream << solution.reconstruct(mesh.globalCellNumber[i], mesh.cells[i]->getCellCenter(), E) << endl;

    wStream << "SCALARS p double" << endl;
    wStream << "LOOKUP_TABLE default" << endl;

    for (int i = 0; i < mesh.nRealCells; ++i)
        wStream << physics.getPressure(solution.reconstruct(mesh.globalCellNumber[i], mesh.cells[i]->getCellCenter())) << endl;


    wStream << "VECTORS U double" << endl;

    for (int i = 0; i < mesh.nRealCells; ++i)
    {
        double rho = solution.reconstruct(mesh.globalCellNumber[i], mesh.cells[i]->getCellCenter(), RHO);
        wStream << solution.reconstruct(mesh.globalCellNumber[i], mesh.cells[i]->getCellCenter(), RHOU) / rho << endl;
        wStream << solution.reconstruct(mesh.globalCellNumber[i], mesh.cells[i]->getCellCenter(), RHOV) / rho << endl;
        wStream << 0.0 << endl;
    }

    //indicator.writeTroubledCellsVTK();

    // get point data

    wStream << "POINT_DATA " << polySize << endl;

    wStream << "SCALARS rho double" << endl;
    wStream << "LOOKUP_TABLE default" << endl;

    for (int i = 0; i < mesh.nRealCells; ++i)
        for (int j = 0; j < mesh.cells[i]->nEntities; ++j)
            wStream << solution.reconstruct(mesh.globalCellNumber[i], *(mesh.cells[i]->nodes[j]), RHO) << endl;

    wStream << "SCALARS e double" << endl;
    wStream << "LOOKUP_TABLE default" << endl;

    for (int i = 0; i < mesh.nRealCells; ++i)
        for (int j = 0; j < mesh.cells[i]->nEntities; ++j)
            wStream << solution.reconstruct(mesh.globalCellNumber[i], *(mesh.cells[i]->nodes[j]), E) << endl;

    wStream << "SCALARS p double" << endl;
    wStream << "LOOKUP_TABLE default" << endl;

    for (int i = 0; i < mesh.nRealCells; ++i)
        for (int j = 0; j < mesh.cells[i]->nEntities; ++j)
            wStream << physics.getPressure(solution.reconstruct(mesh.globalCellNumber[i], *(mesh.cells[i]->nodes[j]))) << endl;


    wStream << "VECTORS U double" << endl;

    for (int i = 0; i < mesh.nRealCells; ++i)
        for (int j = 0; j < mesh.cells[i]->nEntities; ++j)
        {
            double rho = solution.reconstruct(mesh.globalCellNumber[i], *(mesh.cells[i]->nodes[j]), RHO);
            wStream << solution.reconstruct(mesh.globalCellNumber[i], *(mesh.cells[i]->nodes[j]), RHOU) / rho << endl;
            wStream << solution.reconstruct(mesh.globalCellNumber[i], *(mesh.cells[i]->nodes[j]), RHOV) / rho << endl;
            wStream << 0.0 << endl;
        }
}