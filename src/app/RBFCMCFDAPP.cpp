#include "MQBasis.hpp"
#include "controlData.hpp"
#include "domainData.hpp"
#include "json.h"
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
    // build control data
    // ****************************************************************************
    std::shared_ptr<controlData> myControlData;
    if (result.count("case"))
        myControlData =
            std::make_shared<controlData>(result["case"].as<std::string>());
    else
        myControlData = std::make_shared<controlData>();

    // ****************************************************************************
    // Build domain data
    // ****************************************************************************
    auto myDomainData = std::make_shared<DomainData>(myControlData);

    // ****************************************************************************
    // Define RBF Type
    // ****************************************************************************
    auto myRBFBasis = std::make_shared<MQBasis>(myControlData);

    // ****************************************************************************
    // Define solver
    // ****************************************************************************
    SimulationFlow mySimulationFlow(myControlData, myDomainData, myRBFBasis);
    mySimulationFlow.solveDomain();

    std::cout << "test ok" << std::endl;
    return 0;
}