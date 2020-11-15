#ifndef READFRAOMMSH_HPP
#define READFRAOMMSH_HPP
// https://www.manpagez.com/info/gmsh/gmsh-2.8.4/gmsh_56.php
// https://blog.csdn.net/hjq376247328/article/details/47040133
// http://blog.sciencenet.cn/home.php?mod=space&uid=441611&do=blog&quickforward=1&id=758751
#include "dataStructure.hpp"
#include "enumMap.hpp"

#include <Eigen/Dense>
#include <fstream>
#include <string>

bool readFromMsh(const std::string& filePath, std::vector<vec3d<double>>& nodes,
                 std::vector<vec3d<double>>& normals,
                 std::map<std::string, std::vector<size_t>>& groupToNodesMap);

void parseNodes(std::ifstream& fileStream, std::vector<vec3d<double>>& nodes,
                std::map<size_t, size_t>& nodesIndexMap);

void parseElements(std::ifstream& fileStream,
                   std::vector<std::vector<size_t>>& elementNodes,
                   std::vector<elementType>& elementTypes,
                   std::map<size_t, std::vector<size_t>>& groupIndexToNodesMap);

void parseGroupNames(std::ifstream& fileStream,
                     std::map<size_t, std::string>& groupIndexToNameMap);

inline void passWhiteSpace(std::ifstream& fileStream)
{
    char next = fileStream.peek();
    while (next == '\n' || next == ' ' || next == '\t' || next == '\r')
    {
        fileStream.get();
        next = fileStream.peek();
    }
};

inline size_t getNodesNumOnElement(size_t elementTypeTag)
{
    size_t nodesNumOnElement = 0;
    switch (elementTypeTag)
    {
        case 1:
            nodesNumOnElement = 2;  // line
            break;
        case 2:
            nodesNumOnElement = 3;  // TRIANGLE
            break;
        case 3:
            nodesNumOnElement = 4;  // Quad
            break;
        case 4:
            nodesNumOnElement = 4;  // Tetras
            break;
        case 5:
            nodesNumOnElement = 8;  // hexahedron
            break;
        default:
            nodesNumOnElement = -1;
            break;
    }
    return nodesNumOnElement;
};

inline elementType getElementType(size_t elementTypeTag)
{
    elementType type;
    switch (elementTypeTag)
    {
        case 1:
            type = elementType::LINE;
            break;
        case 2:
            type = elementType::TRIANGLE;

            break;
        case 3:
            type = elementType::QUADRANGLE;

            break;
        case 4:
            type = elementType::TETRAHEDRON;
            break;
        case 5:
            type = elementType::HEXAHEDRON;
            break;
        default:
            type = elementType::NONE;
            break;
    }
    return type;
};

#endif
