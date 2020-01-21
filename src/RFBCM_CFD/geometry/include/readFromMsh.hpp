#ifndef READFRAOMMSH_HPP
#define READFRAOMMSH_HPP
// https://www.manpagez.com/info/gmsh/gmsh-2.8.4/gmsh_56.php
// https://blog.csdn.net/hjq376247328/article/details/47040133
// http://blog.sciencenet.cn/home.php?mod=space&uid=441611&do=blog&quickforward=1&id=758751
#include "dataStructure.hpp"
#include <Eigen/Dense>
#include <fstream>
#include <string>

bool readFromMsh(const std::string& filePath, std::vector<vec3d<double>>& nodes,
                 std::map<std::string, std::vector<int>>& groupToNodesMap);

void parseNodes(std::ifstream& fileStream, std::vector<vec3d<double>>& nodes,
                std::map<int, int>& nodesMap);

void parseElements(std::ifstream& fileStream,
                   std::map<int, std::vector<int>>& groupIndexToNodesMap);

void parsePhysicalNames(std::ifstream& fileStream,
                        std::map<int, std::string>& physicalIndexToNameMap);

inline void passWhiteSpace(std::ifstream& fileStream)
{
    char next = fileStream.peek();
    while (next == '\n' || next == ' ' || next == '\t' || next == '\r')
    {
        fileStream.get();
        next = fileStream.peek();
    }
};

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
