#ifndef READFRAOMMSH_HPP
#define READFRAOMMSH_HPP

#include "dataStructure.hpp"
#include <Eigen/Dense>
#include <string>

bool readFromMsh(const std::string& filePath,
                 std::vector<vec3d<double>>& nodes);

void parseNodes(std::ifstream& fileStream, std::vector<vec3d<double>>& nodes);

inline size_t getNumNodesPerElement(int elementType)
{
    int numNodesPerElement = 0;
    switch (elementType)
    {
        case 2:
            numNodesPerElement = 3;  // Triangle
            break;
        case 3:
            numNodesPerElement = 4;  // Quad
            break;
        case 4:
            numNodesPerElement = 4;  // Tet
            break;
        case 5:
            numNodesPerElement = 8;  // hexahedron
            break;
        default:
            numNodesPerElement = -1;
            break;
    }
    return numNodesPerElement;
};

#endif
