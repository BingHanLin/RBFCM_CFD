#include "json.h"
#include "rectangle.hpp"
#include <fstream>
#include <iostream>

int main() {
    // ****************************************************************************
    // Read Json file
    // ****************************************************************************
    std::ifstream  json_stream( "params.json" );
    nlohmann::json myParams;
    json_stream >> myParams;

    // ****************************************************************************
    // Read mesh from inp file
    // ****************************************************************************
    // Create Mesh
    RECTANGLE myMesh( myParams[ "GeometryControl" ][ "NodeNum" ][ 0 ],
                      myParams[ "GeometryControl" ][ "NodeNum" ][ 1 ],
                      myParams[ "GeometryControl" ][ "Size" ][ 0 ],
                      myParams[ "GeometryControl" ][ "Size" ][ 1 ] );

    // Chose Basis
    // GMMMQBasis2D RBFBasis( 6.0 );
    // ****************************************************************************
    // Define Solver
    // ****************************************************************************

    std::cout << "test ok" << std::endl;
    return 0;
}