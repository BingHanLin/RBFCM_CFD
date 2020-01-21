#include "readFromMsh.hpp"
#include "messages.hpp"
#include <iostream>
#include <map>
#include <sstream>
#include <vector>
// https://www.manpagez.com/info/gmsh/gmsh-2.8.4/gmsh_56.php
// http://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format-version-2-_0028Legacy_0029
bool readFromMsh(const std::string& filePath, std::vector<vec3d<double>>& nodes,
                 std::map<std::string, std::vector<int>>& groupToNodesMap)
{
    const auto invalidFormat = [](std::string messages) -> bool {
        ASSERT("Invalid format: " << messages);
        return false;
    };
    const auto notImplemented = [](std::string messages) -> bool {
        ASSERT("Not implemented:" << messages);
        return false;
    };

    // =========================================================================
    // open file and check
    // =========================================================================
    nodes.clear();
    groupToNodesMap.clear();

    std::ifstream fileStream(filePath.c_str(), std::ios::in);

    if (!fileStream.is_open())
    {
        ASSERT("msh file is not found in path: " << filePath);
        return false;
    }

    // =========================================================================
    // parse header and check
    // =========================================================================
    std::string buf;
    fileStream >> buf;
    if (buf != "$MeshFormat")
    {
        return invalidFormat("msh MeshFormat tag not found");
    }

    double versionNumber;
    int fileType;  // equal to 0 indicates the ASCII file format.
    size_t dataSize;
    fileStream >> versionNumber >> fileType >> dataSize;

    if (fileType != 0)
    {
        return notImplemented("only accept ASCII file format for .msh");
    }

    if (dataSize != 8)
    {
        return invalidFormat("data size must be 8 bytes");
    }

    fileStream >> buf;
    if (buf != "$EndMeshFormat")
    {
        return invalidFormat("msh EndMeshFormat tag not found");
    }

    // =========================================================================
    // parse other tags
    // =========================================================================
    std::map<int, int> nodesMap;
    std::map<int, std::vector<int>> groupIndexToNodesMap;
    std::map<int, std::string> groupIndexToNamesMap;

    while (!fileStream.eof())
    {
        buf.clear();
        fileStream >> buf;
        // parse Nodes
        if (buf == "$Nodes")
        {
            parseNodes(fileStream, nodes, nodesMap);
            fileStream >> buf;
            if (buf != "$EndNodes")
            {
                return invalidFormat("msh EndNodes tag not found");
            }
        }
        // parse Elements
        else if (buf == "$Elements")
        {
            parseElements(fileStream, groupIndexToNodesMap);
            fileStream >> buf;
            if (buf != "$EndElements")
            {
                return invalidFormat("msh EndElements tag not found");
            }
            std::cout << "test" << std::endl;
        }
        // parse physicalNames
        else if (buf == "$PhysicalNames")
        {
            parsePhysicalNames(fileStream, groupIndexToNamesMap);
            fileStream >> buf;
            if (buf != "$EndPhysicalNames")
            {
                return invalidFormat("msh EndPhysicalNames tag not found");
            }
            std::cout << "test" << std::endl;
        }
    }
    fileStream.close();

    // =========================================================================
    // processing in mesh data
    // =========================================================================

    //  compact nodes in groupIndexToNodesMap
    for (auto& oneMap : groupIndexToNodesMap)
    {
        std::sort(oneMap.second.begin(), oneMap.second.end());
        oneMap.second.erase(unique(oneMap.second.begin(), oneMap.second.end()),
                            oneMap.second.end());
    }

    //  map node index to std::vector container index
    for (auto& oneMap : groupIndexToNodesMap)
    {
        for (auto& oneNodeIndex : oneMap.second)
        {
            oneNodeIndex = nodesMap.at(oneNodeIndex);
        }
    }

    for (auto& oneMap : groupIndexToNamesMap)
    {
        groupToNodesMap.insert(std::pair<std::string, std::vector<int>>(
            oneMap.second, groupIndexToNodesMap.at(oneMap.first)));
    }

    return true;
}

void parseNodes(std::ifstream& fileStream, std::vector<vec3d<double>>& nodes,
                std::map<int, int>& nodesMap)
{
    size_t numOfNodes;
    fileStream >> numOfNodes;
    nodes.resize(numOfNodes);

    for (size_t i = 0; i < numOfNodes; i++)
    {
        int nodeIndex;
        fileStream >> nodeIndex;

        nodesMap.insert(std::pair<int, int>(nodeIndex, i));

        fileStream >> nodes[i](0) >> nodes[i](1) >> nodes[i](2);
    }
};

void parseElements(std::ifstream& fileStream,
                   std::map<int, std::vector<int>>& groupIndexToNodesMap)
{
    size_t numOfElements;
    fileStream >> numOfElements;

    for (size_t i = 0; i < numOfElements; i++)
    {
        // parse per element header
        int elementIndex, elementType, numTags;
        fileStream >> elementIndex >> elementType >> numTags;

        // https://www.manpagez.com/info/gmsh/gmsh-2.8.4/gmsh_56.php
        int groupIndex;
        for (size_t j = 0; j < numTags; j++)
        {
            int tag;
            fileStream >> tag;

            if (j == 0)
            {
                groupIndex = tag;
            }
        }

        // parse node index per element.
        std::vector<int> nodesOfElements;
        size_t numNodesPerElement = getNumNodesPerElement(elementType);
        for (size_t j = 0; j < numNodesPerElement; j++)
        {
            int nodeIndex;
            fileStream >> nodeIndex;
            nodesOfElements.push_back(nodeIndex);
        }

        if (groupIndexToNodesMap.find(groupIndex) == groupIndexToNodesMap.end())
        {
            groupIndexToNodesMap.insert(
                std::pair<int, std::vector<int>>(groupIndex, nodesOfElements));
        }
        else
        {
            groupIndexToNodesMap.at(groupIndex)
                .insert(groupIndexToNodesMap.at(groupIndex).end(),
                        nodesOfElements.begin(), nodesOfElements.end());
        }
    }
};

void parsePhysicalNames(std::ifstream& fileStream,
                        std::map<int, std::string>& groupIndexToNamesMap)
{
    size_t numOfGroups;
    fileStream >> numOfGroups;
    for (size_t i = 0; i < numOfGroups; i++)
    {
        int dimension, groupIndex;
        std::string groupName;
        fileStream >> dimension >> groupIndex >> groupName;
        groupIndexToNamesMap.insert(
            std::pair<int, std::string>(groupIndex, groupName));
    }
};