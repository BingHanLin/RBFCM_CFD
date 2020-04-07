#include "meshData.hpp"
#include "constantValueBC.hpp"
#include "readFromMsh.hpp"
#include "rectangle.hpp"

#include <iostream>

MeshData::MeshData(nlohmann::json& controlParams)
    : geometryControlParams_(controlParams.at("GeometryControl")),
      physicsControlParams_(controlParams.at("PhysicsControl")),
      numOfNodes_(0),
      groupToNodesMap_(),
      nodes_(),
      normals_()
{
    std::map<std::string, std::vector<int>> groupToNodesMapBeforeCompact;

    if (geometryControlParams_.at("Type") == meshTypeEnum::DEFAULT)
    {
        std::cout << "read nodes from msh file" << std::endl;
        bool isReadSuccess =
            readFromMsh(geometryControlParams_.at("Path"), nodes_, normals_,
                        groupToNodesMapBeforeCompact);
    }
    else if (geometryControlParams_.at("Type") == meshTypeEnum::RECTNAGLE)
    {
        std::cout << "create rectangle mesh and read from msh file"
                  << std::endl;
    }
    else
    {
        std::cout << "mesh type is not defined" << std::endl;
    }

    numOfNodes_ = nodes_.size();

    // compactGroupToNodesMap();
    buildBoundaryConditions();
    matchBoundaryConditions(groupToNodesMapBeforeCompact);
}

void MeshData::buildBoundaryConditions()
{
    auto& constantValueBCData =
        physicsControlParams_.at("Boundary").at("ConstantValue");
    for (auto& oneConstantValueBCData : constantValueBCData)
    {
        const std::string groupName = oneConstantValueBCData.at("groupName");
        auto constantValueBC = std::make_shared<ConstantValueBC>(
            oneConstantValueBCData.at("value"));

        groupToBCMap_.insert({groupName, constantValueBC});
    }
}

void MeshData::matchBoundaryConditions(
    std::map<std::string, std::vector<int>> groupToNodesMapBeforeCompact)
{
    nodesBC_.resize(nodes_.size());
    auto& constantValueBCData =
        physicsControlParams_.at("Boundary").at("ConstantValue");
    for (auto& oneConstantValueBCData : constantValueBCData)
    {
        const std::string groupNmae = oneConstantValueBCData.at("groupName");
        for (int& nodeIndex : groupToNodesMapBeforeCompact.at(groupNmae))
        {
            nodesBC_[nodeIndex] = groupToBCMap_.at(groupNmae);
        }
    }
}

int MeshData::numOfNodes() const
{
    return numOfNodes_;
};
