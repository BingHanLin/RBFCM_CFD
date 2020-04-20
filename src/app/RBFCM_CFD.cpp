#include "MQBasis.hpp"
#include "controlData.hpp"
#include "json.h"
#include "meshData.hpp"
#include "simulationDomain.hpp"

#include <fstream>
#include <iostream>

int main()
{
    // ****************************************************************************
    // build control data
    // ****************************************************************************
    controlData::instance();

    // ****************************************************************************
    // Build mesh data
    // ****************************************************************************
    auto myMeshData = std::make_shared<MeshData>();

    // ****************************************************************************
    // Define RBF Type
    // ****************************************************************************
    auto myRBFBasis = std::make_shared<MQBasis>();

    // ****************************************************************************
    // Define Solver
    // ****************************************************************************
    SimulationDomain mySimulationDomain(myMeshData, myRBFBasis);
    // mySimulationDomain.solveDomain();

    // mySimulationDomain.exportData();

    std::cout << "test ok" << std::endl;
    return 0;
}