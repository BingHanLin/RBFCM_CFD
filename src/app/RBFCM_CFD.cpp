#include "MQBasis.hpp"
#include "json.h"
#include "meshData.hpp"
#include "simulationDomain.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>

int main()
{
    // ****************************************************************************
    // Read Json file
    // ****************************************************************************
    std::ifstream json_stream("params.json");
    nlohmann::json myParams;
    json_stream >> myParams;

    // ****************************************************************************
    // Build mesh data
    // ****************************************************************************
    auto myMeshData = std::make_shared<MeshData>(myParams);

    // ****************************************************************************
    // Define RBF Type
    // ****************************************************************************
    auto myRBFBasis = std::make_shared<MQBasis>(
        myParams.at("SolverConstrol").at("RBFCoefficient"));

    // ****************************************************************************
    // Define Solver
    // ****************************************************************************
    SimulationDomain mySimulationDomain(myMeshData, myRBFBasis, myParams);
    mySimulationDomain.solveDomain();
    // mySimulationDomain.exportData();

    std::cout << "test ok" << std::endl;
    return 0;
}