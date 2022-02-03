
#include "KDTreeEigenAdaptor.hpp"
#include "catch.hpp"
#include "dataStructure.hpp"
#include <iostream>

TEST_CASE("KDTreeEigenAdaptor Test", "[RBF]")
{
    SECTION("KDTreeEigenAdaptor")
    {
        std::vector<vec3d<double>> nodes = {
            {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0},
            {0.0, 2.0, 0.0}, {0.0, 0.0, 3.0}, {0.0, 5.0, 0.0},
            {3.0, 4.0, 0.0}, {0.0, 5.1, 0.0}, {5.0, 5.0, 0.0}};

        const double neighborRadius = 5.0;

        const auto kdTree =
            KDTreeEigenAdaptor<std::vector<vec3d<double>>, double, 3>(nodes);

        std::vector<size_t> neighboursID;
        kdTree.query(0, neighborRadius, neighboursID);

        REQUIRE(neighboursID.size() == 50);
    }
}