#include "readFromMsh.hpp"
#include "messages.hpp"
#include <iostream>
#include <map>
#include <sstream>
#include <vector>
// https://www.manpagez.com/info/gmsh/gmsh-2.8.4/gmsh_56.php
// http://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format-version-2-_0028Legacy_0029
bool readFromMsh(const std::string& filePath, std::vector<vec3d<double>>& nodes,
                 std::vector<vec3d<double>>& normals,
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
    std::map<int, int> nodesIndexMap;
    std::vector<std::vector<int>> elementNodes;
    std::vector<int> elementTypes;
    std::map<int, std::vector<int>> groupIndexToNodesMap;
    std::map<int, std::string> groupIndexToNamesMap;

    while (!fileStream.eof())
    {
        buf.clear();
        fileStream >> buf;
        // parse Nodes
        if (buf == "$Nodes")
        {
            std::cout << "parse Nodes" << std::endl;
            parseNodes(fileStream, nodes, nodesIndexMap);
            fileStream >> buf;
            if (buf != "$EndNodes")
            {
                return invalidFormat("msh EndNodes tag not found");
            }
        }
        // parse Elements
        else if (buf == "$Elements")
        {
            std::cout << "parse Elements" << std::endl;
            parseElements(fileStream, elementNodes, elementTypes,
                          groupIndexToNodesMap);

            fileStream >> buf;
            if (buf != "$EndElements")
            {
                return invalidFormat("msh EndElements tag not found");
            }
        }
        // parse physicalNames
        else if (buf == "$PhysicalNames")
        {
            std::cout << "parse PhysicalNames" << std::endl;
            parseGroupNames(fileStream, groupIndexToNamesMap);
            fileStream >> buf;
            if (buf != "$EndPhysicalNames")
            {
                return invalidFormat("msh EndPhysicalNames tag not found");
            }
        }
    }
    fileStream.close();

    // =========================================================================
    // processing in mesh data
    // =========================================================================
    std::cout << "parse PhysicalNamesprocessing in mesh data..." << std::endl;

    //  compact nodes in groupIndexToNodesMap
    std::cout << "    compact nodes in groupIndexToNodesMap" << std::endl;
    for (auto& oneMap : groupIndexToNodesMap)
    {
        std::sort(oneMap.second.begin(), oneMap.second.end());
        oneMap.second.erase(unique(oneMap.second.begin(), oneMap.second.end()),
                            oneMap.second.end());
    }

    //  compute normals on nodes by loop elements
    normals.resize(nodes.size());
    for (int i = 0; i < elementNodes.size(); i++)
    {
        if (elementTypes[i] == elementType::LINE)
        {
            int n0 = nodesIndexMap.at(elementNodes[i][0]);
            int n1 = nodesIndexMap.at(elementNodes[i][1]);

            vec3d<double> v0 = nodes[n0];
            vec3d<double> v1 = nodes[n1];

            vec3d<double> normal(-(v1(1) - v0(1)), v1(0) - v0(0), 0);
            normal.normalize();

            normals[n0] += normal;
            normals[n1] += normal;
        }
        // if (oneElementNodes.size() == 3)
        // {
        //     int n0 = nodesIndexMap.at(oneElementNodes[0]);
        //     int n1 = nodesIndexMap.at(oneElementNodes[1]);
        //     int n2 = nodesIndexMap.at(oneElementNodes[2]);

        //     vec3d<double> v0 = nodes[n0];
        //     vec3d<double> v1 = nodes[n1];
        //     vec3d<double> v2 = nodes[n2];

        //     vec3d<double> normal = (v1 - v0).cross(v2 - v0);
        //     normal.normalize();

        //     normals[n0] += normal;
        //     normals[n1] += normal;
        //     normals[n2] += normal;
        // }
        // else if (oneElementNodes.size() == 4)
        // {
        //     int n0 = nodesIndexMap.at(oneElementNodes[0]);
        //     int n1 = nodesIndexMap.at(oneElementNodes[1]);
        //     int n2 = nodesIndexMap.at(oneElementNodes[2]);
        //     int n3 = nodesIndexMap.at(oneElementNodes[3]);

        //     vec3d<double> v0 = nodes[n0];
        //     vec3d<double> v1 = nodes[n1];
        //     vec3d<double> v2 = nodes[n2];
        //     vec3d<double> v3 = nodes[n3];

        //     vec3d<double> normal = (v1 - v0).cross(v3 - v0);
        //     normal.normalize();

        //     normals[n0] += normal;
        //     normals[n1] += normal;
        //     normals[n2] += normal;
        //     normals[n3] += normal;
        // }
    }

    // normalize all normals
    std::for_each(normals.begin(), normals.end(),
                  [](vec3d<double>& normal) { normal.normalize(); });

    // std::cout << "node size: " << nodes.size() << std::endl;
    // std::for_each(normals.begin(), normals.end(), [](vec3d<double>& normal) {
    //     std::cout << normal(0) << "; " << normal(1) << "; " << normal(2) <<
    //     "; "
    //               << std::endl;
    // });

    //  map node index to std::vector container index
    std::cout << "    map node index to std::vector container index"
              << std::endl;
    for (auto& oneMap : groupIndexToNodesMap)
    {
        for (auto& oneNodeIndex : oneMap.second)
        {
            oneNodeIndex = nodesIndexMap.at(oneNodeIndex);
        }
    }

    //  map node index to groupToNodesMap
    std::cout << "    map node index to groupToNodesMap" << std::endl;
    for (auto& oneMap : groupIndexToNamesMap)
    {
        groupToNodesMap.insert(std::pair<std::string, std::vector<int>>(
            oneMap.second, groupIndexToNodesMap.at(oneMap.first)));
    }

    return true;
}

void parseNodes(std::ifstream& fileStream, std::vector<vec3d<double>>& nodes,
                std::map<int, int>& nodesIndexMap)
{
    size_t numOfNodes;
    fileStream >> numOfNodes;
    nodes.resize(numOfNodes);

    for (size_t i = 0; i < numOfNodes; i++)
    {
        int nodeIndex;
        fileStream >> nodeIndex;
        nodesIndexMap.insert(std::pair<int, int>(nodeIndex, i));
        fileStream >> nodes[i](0) >> nodes[i](1) >> nodes[i](2);
    }
};

void parseElements(std::ifstream& fileStream,
                   std::vector<std::vector<int>>& elementNodes,
                   std::vector<int>& elementTypes,
                   std::map<int, std::vector<int>>& groupIndexToNodesMap)
{
    size_t numOfElements;
    fileStream >> numOfElements;

    elementNodes.resize(numOfElements);
    elementTypes.resize(numOfElements);

    for (size_t i = 0; i < numOfElements; i++)
    {
        // parse per element header
        int elementIndex, elementType, numTags;
        fileStream >> elementIndex >> elementType >> numTags;

        elementTypes[i] = getElementType(elementType);

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
        size_t nodesNumOnElement = getNodesNumOnElement(elementType);
        std::vector<int> nodesOnElement(nodesNumOnElement);
        for (size_t j = 0; j < nodesNumOnElement; j++)
        {
            int nodeIndex;
            fileStream >> nodeIndex;
            nodesOnElement[j] = nodeIndex;
        }
        elementNodes[i] = nodesOnElement;

        if (groupIndexToNodesMap.find(groupIndex) == groupIndexToNodesMap.end())
        {
            groupIndexToNodesMap.insert(
                std::pair<int, std::vector<int>>(groupIndex, nodesOnElement));
        }
        else
        {
            groupIndexToNodesMap.at(groupIndex)
                .insert(groupIndexToNodesMap.at(groupIndex).end(),
                        nodesOnElement.begin(), nodesOnElement.end());
        }
    }
};

void parseGroupNames(std::ifstream& fileStream,
                     std::map<int, std::string>& groupIndexToNamesMap)
{
    size_t numOfGroups;
    fileStream >> numOfGroups;
    for (size_t i = 0; i < numOfGroups; i++)
    {
        int dimension, groupIndex;
        std::string groupNameWithQuotes;
        fileStream >> dimension >> groupIndex >> groupNameWithQuotes;

        std::istringstream groupNameStream(groupNameWithQuotes);
        std::string skip;
        std::string groupName;
        std::getline(std::getline(groupNameStream, skip, '"'), groupName, '"');

        groupIndexToNamesMap.insert(
            std::pair<int, std::string>(groupIndex, groupName));
    }
};