#include "MQBasis.hpp"
#include "MeshData.hpp"
#include "controlData.hpp"
#include "scalarConditionPool.hpp"
#include "simulationFlow.hpp"

#include <cxxopts.hpp>
#include <fstream>
#include <iostream>

int main(int argc, char** argv)
{
    // ****************************************************************************
    // parse arguemnt
    // ****************************************************************************
    cxxopts::Options options("RBFCM CFD",
                             "Radial basis function collocation method (RBFCM) "
                             "for computational fluid dynamics");

    options.add_options()("c,case", "Case Path", cxxopts::value<std::string>())(
        "h,help", "Print usage");

    auto result = options.parse(argc, argv);

    if (result.count("help"))
    {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    // ****************************************************************************
    // Build control data
    // ****************************************************************************
    std::shared_ptr<ControlData> myControlData;
    if (result.count("case"))
        myControlData =
            std::make_shared<ControlData>(result["case"].as<std::string>());
    else
        myControlData = std::make_shared<ControlData>();

    // ****************************************************************************
    // Build mesh data
    // ****************************************************************************
    auto myMeshData = std::make_shared<MeshData>(myControlData.get());

    // ****************************************************************************
    // Define RBF type
    // ****************************************************************************
    auto myRBFBasis =
        std::make_shared<MQBasis>(myControlData.get(), myMeshData.get());

    // ****************************************************************************
    // Build condition pool
    // ****************************************************************************
    auto myConditionPoolPool = std::make_shared<ScalarConditionPool>(
        myControlData.get(), myMeshData.get());

    // ****************************************************************************
    // Define solver
    // ****************************************************************************
    SimulationFlow mySimulationFlow(myControlData.get(), myRBFBasis.get(),
                                    myMeshData.get(),
                                    myConditionPoolPool.get());
    mySimulationFlow.solveDomain();

    std::cout << "test ok" << std::endl;
    return 0;
}