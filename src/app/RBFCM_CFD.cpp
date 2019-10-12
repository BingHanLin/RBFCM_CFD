#include "MQBasis2D.hpp"
#include "json.h"
#include "jsonEnumMap.hpp"
#include "rectangle.hpp"
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
    // Read mesh from inp file
    // ****************************************************************************
    // nlohmann::json meshType = myParams[ "GeometryControl" ][ "Type" ];
    // meshTypeEnum   test     = meshType.get< meshTypeEnum >();

    typedef Rectangle myMeshType;
    Rectangle myMesh(myParams["GeometryControl"]["NodeNum"][0],
                     myParams["GeometryControl"]["NodeNum"][1],
                     myParams["GeometryControl"]["Size"][0],
                     myParams["GeometryControl"]["Size"][1]);

    // ****************************************************************************
    // Define RBF Type
    // ****************************************************************************
    typedef MQBasis2D myRBFBasisType;
    MQBasis2D myRBFBasis;

    // ****************************************************************************
    // Define Solver
    // ****************************************************************************
    SimulationDomain<myMeshType, myRBFBasisType> mySimulationDomain(
        myMesh, myRBFBasis, myParams);

    std::cout << myRBFBasis.cc_ << std::endl;
    std::cout << "test ok" << std::endl;
    return 0;
}