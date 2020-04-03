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
    MeshData myMeshData(myParams.at("GeometryControl"));

    // ****************************************************************************
    // Define RBF Type
    // ****************************************************************************
    MQBasis myRBFBasis(myParams["SolverConstrol"]["RBFCoefficient"]);

    // ****************************************************************************
    // Define Solver
    // ****************************************************************************
    // SimulationDomain<myMeshType, myRBFBasisType> mySimulationDomain(
    //     myMesh, myRBFBasis, myParams);

    // mySimulationDomain.solveDomain();
    // mySimulationDomain.exportData();

    std::cout << "test ok" << std::endl;
    return 0;
}