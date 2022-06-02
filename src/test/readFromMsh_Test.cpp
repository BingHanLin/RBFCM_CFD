#include "catch.hpp"
#include "messages.hpp"
#include "readFromMsh.hpp"
#include <filesystem>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

TEST_CASE("readFromMsh Test", "[readFromMsh]")
{
    SECTION("readFromMsh")
    {
        const std::string filePath = std::filesystem::current_path()
                                         .concat("/data/2dRec10.msh")
                                         .string();

        std::vector<vec3d<double>> nodes;
        std::vector<vec3d<double>> normals;
        std::map<std::string, std::vector<size_t>> groupToNodesMap;

        readFromMsh(filePath, nodes, normals, groupToNodesMap);

        REQUIRE(nodes.size() == 100);
        REQUIRE(nodes.size() == normals.size());
        REQUIRE(groupToNodesMap.find("left") != groupToNodesMap.end());
        REQUIRE(groupToNodesMap.find("top") != groupToNodesMap.end());
        REQUIRE(groupToNodesMap.find("right") != groupToNodesMap.end());
        REQUIRE(groupToNodesMap.find("bottom") != groupToNodesMap.end());
        REQUIRE(groupToNodesMap.find("inner") != groupToNodesMap.end());
    }
}
