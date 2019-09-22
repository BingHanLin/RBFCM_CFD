#include "GMMMQBasis2D.h"
#include "gmm.h"
#include "json.h"
#include <fstream>
#include <iostream>

int main() {
    // ****************************************************************************
    // Read Json file
    // ****************************************************************************
    std::ifstream  json_stream( "params.json" );
    nlohmann::json my_params;
    json_stream >> my_params;

    // ****************************************************************************
    // Read mesh from inp file
    // ****************************************************************************
    // Create Mesh
    GMMRECTANGLE MESH( m, n, Lx, Ly );
    // Chose Basis
    GMMMQBasis2D RBFBasis( 6.0 );
    // ****************************************************************************
    // Define Solver
    // ****************************************************************************

    std::cout << "test ok" << std::endl;
    return 0;
}