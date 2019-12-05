#include "MQBasis2D.hpp"
#include "json.h"
#include "meshData.hpp"
#include "simulationDomain.hpp"
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
    // // Define RBF Type
    // //
    // ****************************************************************************
    // typedef MQBasis2D myRBFBasisType;
    // MQBasis2D myRBFBasis(myParams["SolverConstrol"]["RBFCoefficient"]);

    // //
    // ****************************************************************************
    // // Define Solver
    // //
    // ****************************************************************************
    // SimulationDomain<myMeshType, myRBFBasisType> mySimulationDomain(
    //     myMesh, myRBFBasis, myParams);

    // mySimulationDomain.solveDomain();
    // mySimulationDomain.exportData();

    std::cout << "test ok" << std::endl;
    return 0;
}