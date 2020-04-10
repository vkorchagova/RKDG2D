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
    cout << "Exporting VTK frame " << fileName << "..." << endl;
    
    ofstream wStream (fileName);

    cout << "\tWriting VTK mesh... ";
    exportMeshVTK(wStream);
    cout << "OK" << endl;

    exportFrameVTK(wStream);
    wStream.close();

    cout << "Exporting VTK frame OK" << endl;
    
    cout << endl;
}


void Writer::exportNativeCoeffs(const std::string& fileName) const
{
    ofstream wStream (fileName);
    cout << "Export native coeffs... ";
    outputFullNativeCoeffs(wStream);
    cout << "OK" << endl;
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

    //  for (const shared_ptr<Cell> cell : mesh.cells)
    cout << "\tWriting rho field... ";
    for (int iCell = 0; iCell < mesh.nRealCells; ++iCell)
        wStream << solution.solToExport[iCell][0] << endl;
    cout << "OK" << endl;

    wStream << "SCALARS e double" << endl;
    wStream << "LOOKUP_TABLE default" << endl;

    cout << "\tWriting e field... ";
    for (int iCell = 0; iCell < mesh.nRealCells; ++iCell)
        wStream << solution.solToExport[iCell][4] << endl;
    cout << "OK" << endl;

    wStream << "SCALARS p double" << endl;
    wStream << "LOOKUP_TABLE default" << endl;

    cout << "\tWriting p field... ";
    for (int iCell = 0; iCell < mesh.nRealCells; ++iCell)
        wStream << solution.solToExport[iCell][5] << endl;
    cout << "OK" << endl;


    wStream << "VECTORS U double" << endl;

    cout << "\tWriting U field... ";
    for (int iCell = 0; iCell < mesh.nRealCells; ++iCell)
    {
        wStream << solution.solToExport[iCell][1] << endl;
        wStream << solution.solToExport[iCell][2] << endl;
        wStream << 0.0 << endl;
    }
    cout << "OK" << endl;

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
