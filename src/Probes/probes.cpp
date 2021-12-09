#include "probes.h"

Probes::Probes(const Mesh& msh, const Solution& sln, const Physics& phs, const std::vector<Point>& _probes) : mesh(msh), solution(sln), physics(phs), probes(_probes)
{
    iCellsForProbes.resize(probes.size());
    pForProbes.resize(probes.size());
    pForProbesRecv.resize(probes.size());

    for (int iProbe = 0; iProbe < probes.size(); ++iProbe)
    {
        iCellsForProbes[iProbe] = -1.0;
        pForProbes[iProbe] = -1.0;
        pForProbesRecv[iProbe] = -1.0;

        for (int iCell = 0; iCell < mesh.cells.size(); ++iCell)
        {
            if (mesh.cells[iCell]->insideCell(probes[iProbe]))
            {
                iCellsForProbes[iProbe] = iCell;
                break;
            }
        }
    }

    if (myRank == 0)
    {
        probesStream.open("probes");
    }
}

Probes::~Probes()
{
    if (myRank == 0)
        probesStream.close();
        
}

void Probes::writeProbes(double t)
{
    for (int iProbe = 0; iProbe < probes.size(); ++iProbe)
    {
        pForProbes[iProbe] = -1;

        if (iCellsForProbes[iProbe] > 0)
        {
            // std::cout << "myRank = " << myRank << ", iCell = " << iCellsForProbes[iProbe] << " probes[iProbe] = " << probes[iProbe] <<  std::endl;
            numvector<double, dimPh> solProbe = solution.reconstruct(iCellsForProbes[iProbe], probes[iProbe]);
            //  std::cout << "myRank = " << myRank << ", solProbe = " << solProbe << std::endl;
            pForProbes[iProbe] = physics.getPressure(solProbe);
            // std::cout << "myRank = " << myRank << ", iCell = " << iCellsForProbes[iProbe] << ", p = " << pForProbes[iProbe] << std::endl; 
        }
    }

    // if (myRank == 0 || myRank == 29)
    // {
    // //     // std::ofstream probesStream.open("probes", std::ios_base::app);

    //     // std::cout << "before " << myRank << ' ' << t << ' ';

    //     for (int iProbe = 0; iProbe < probes.size(); ++iProbe)
    //     {
    //         std::cout << pForProbes[iProbe] << " vs " << pForProbesRecv[iProbe] << " | ";    
    //     }

    //     std::cout << std::endl;

    // //     // probesStream.close();
    // }

    // MPI_Barrier(MPI_COMM_WORLD);

    // std::cout << "before reduce" << std::endl;
// for (int iProbe = 0; iProbe < probes.size(); ++iProbe)
    // {
        MPI_Reduce(
            &pForProbes[0],
            &pForProbesRecv[0],
            pForProbesRecv.size(),
            MPI_DOUBLE,
            MPI_MAX,
            0,
            MPI_COMM_WORLD
        );
    // }
    

    if (myRank == 0 )
    {
        // std::ofstream probesStream.open("probes", std::ios_base::app);

        probesStream << t << ' ';

        for (int iProbe = 0; iProbe < pForProbesRecv.size(); ++iProbe)
        {
            probesStream  << pForProbesRecv[iProbe] << ' ';   
        }

        probesStream << std::endl;
    }
}