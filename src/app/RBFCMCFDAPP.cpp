#include "MQBasis.hpp"
#include "PUConditionPool.hpp"
#include "controlData.hpp"
#include "enumMap.hpp"
#include "incompressibleDomain.hpp"
#include "meshData.hpp"
#include "scalarConditionPool.hpp"
#include "scalarTransportDomain.hpp"

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
    const double shapeParameter =
        myControlData->paramsDataAt({"solverControl", "shapeParameter"});

    const size_t dimension =
        myControlData->paramsDataAt({"solverControl", "dimension"});

    auto myRBFBasis =
        std::make_shared<MQBasis>(shapeParameter, dimension, myMeshData.get());

    const auto mySolverType =
        myControlData->paramsDataAt({"solverType"}).get<solverType>();

    if (mySolverType == solverType::SCALARTRANSPORT)
    {
        // ****************************************************************************
        // Build condition pool
        // ****************************************************************************
        auto myConditionPoolPool = std::make_shared<ScalarConditionPool>(
            myControlData.get(), myMeshData.get());

        // ****************************************************************************
        // Define solver
        // ****************************************************************************
        ScalarTransportDomain mySolverDomain(myControlData.get(),
                                             myRBFBasis.get(), myMeshData.get(),
                                             myConditionPoolPool.get());
        mySolverDomain.solveDomain();
    }
    else
    {
        // ****************************************************************************
        // Build condition pool
        // ****************************************************************************
        auto myConditionPoolPool = std::make_shared<PUConditionPool>(
            myControlData.get(), myMeshData.get());

        // ****************************************************************************
        // Define solver
        // ****************************************************************************
        IncompressibleDomain mySolverDomain(myControlData.get(),
                                            myRBFBasis.get(), myMeshData.get(),
                                            myConditionPoolPool.get());
        mySolverDomain.solveDomain();
    }

    std::cout << "test ok" << std::endl;
    return 0;
}