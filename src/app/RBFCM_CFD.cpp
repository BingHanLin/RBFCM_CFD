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
    auto myControlData = std::make_shared<controlData>();

    // ****************************************************************************
    // Build mesh data
    // ****************************************************************************
    auto myMeshData = std::make_shared<MeshData>(myControlData);

    // ****************************************************************************
    // Define RBF Type
    // ****************************************************************************
    auto myRBFBasis = std::make_shared<MQBasis>(myControlData);

    // ****************************************************************************
    // Define Solver
    // ****************************************************************************
    SimulationDomain mySimulationDomain(myControlData, myMeshData, myRBFBasis);
    mySimulationDomain.solveDomain();
    // mySimulationDomain.exportData();

    std::cout << "test ok" << std::endl;
    return 0;
}