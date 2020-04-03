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
                 std::map<std::string, std::vector<int>>& groupToNodesMap);

void parseNodes(std::ifstream& fileStream, std::vector<vec3d<double>>& nodes,
                std::map<int, int>& nodesIndexMap);

void parseElements(std::ifstream& fileStream,
                   std::vector<std::vector<int>>& elementNodes,
                   std::vector<int>& elementTypes,
                   std::map<int, std::vector<int>>& groupIndexToNodesMap);

void parseGroupNames(std::ifstream& fileStream,
                     std::map<int, std::string>& groupIndexToNameMap);

inline void passWhiteSpace(std::ifstream& fileStream)
{
    char next = fileStream.peek();
    while (next == '\n' || next == ' ' || next == '\t' || next == '\r')
    {
        fileStream.get();
        next = fileStream.peek();
    }
};

inline size_t getNodesNumOnElement(int elementType)
{
    int nodesNumOnElement = 0;
    switch (elementType)
    {
        case 1:
            nodesNumOnElement = 2;  // line
            break;
        case 2:
            nodesNumOnElement = 3;  // Triangle
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

inline elementType getElementType(int elementTypeTag)
{
    elementType type;
    switch (elementTypeTag)
    {
        case 1:
            type = elementType::Line;
            break;
        case 2:
            type = elementType::Triangle;

            break;
        case 3:
            type = elementType::Quadrangle;

            break;
        case 4:
            type = elementType::Tetrahedron;
            break;
        case 5:
            type = elementType::Hexahedron;
            break;
        default:
            type = elementType::None;
            break;
    }
    return type;
};

#endif
