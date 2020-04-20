#include "meshData.hpp"
#include "constantValueBC.hpp"
#include "readFromMsh.hpp"
#include "rectangle.hpp"

#include <iostream>

MeshData::MeshData()
    : controlData_(controlData::instance()),
      numOfNodes_(),
      groupToNodesMap_(),
      nodes_(),
      normals_(),
      nodesToGroup_()
{
    meshTypeEnum meshType =
        controlData_->paramsDataAt({"GeometryControl", "Type"});

    if (meshType == meshTypeEnum::DEFAULT)
    {
        std::cout << "read nodes from msh file" << std::endl;
        std::map<std::string, std::vector<int>> groupToNodesMapBeforeCompact;

        const std::string meshFileName =
            controlData_->paramsDataAt({"GeometryControl", "FileName"})
                .get<std::string>();

        const std::string absPath =
            controlData_->workingDir().concat("/" + meshFileName).string();

        bool isReadSuccess = readFromMsh(absPath, nodes_, normals_,
                                         groupToNodesMapBeforeCompact);

        compactGroupToNodesMap(groupToNodesMapBeforeCompact);
    }
    else if (meshType == meshTypeEnum::RECTNAGLE)
    {
        std::cout << "create rectangle mesh and read from msh file"
                  << std::endl;
    }
    else
    {
        std::cout << "mesh type is not defined" << std::endl;
    }

    numOfNodes_ = nodes_.size();
    kdTree_ = KDTreeEigenAdaptor<std::vector<vec3d<double>>, double, 3>(nodes_);
    buildBoundaryConditions();
}

void MeshData::compactGroupToNodesMap(
    std::map<std::string, std::vector<int>> groupToNodesMapBeforeCompact)
{
    nodesToGroup_.resize(nodes_.size());
    groupToNodesMap_.clear();

    const auto& constantValueBCData = controlData_->paramsDataAt(
        {"PhysicsControl", "Boundary", "ConstantValue"});

    for (auto& oneBCData : constantValueBCData)
    {
        const std::string groupName = oneBCData.at("GroupName");
        groupToNodesMap_.insert({groupName, {}});

        // TODO: parallel
        for (int& nodeIndex : groupToNodesMapBeforeCompact.at(groupName))
        {
            nodesToGroup_[nodeIndex] = groupName;
        }
    }

    for (int nodeID = 0; nodeID < nodes_.size(); nodeID++)
    {
        groupToNodesMap_[nodesToGroup_[nodeID]].push_back(nodeID);
    }
}

void MeshData::buildBoundaryConditions()
{
    nodesToBC_.resize(nodes_.size());

    const auto& constantValueBCData = controlData_->paramsDataAt(
        {"PhysicsControl", "Boundary", "ConstantValue"});

    for (auto& oneBCData : constantValueBCData)
    {
        const std::string groupName = oneBCData.at("GroupName");

        auto constantValueBC =
            std::make_shared<ConstantValueBC>(oneBCData.at("Value"));
        groupToBCMap_.insert({groupName, constantValueBC});

        for (int& nodeID : groupToNodesMap_.at(groupName))
        {
            nodesToBC_[nodeID] = constantValueBC;
        }
    }
}

nodesCloud MeshData::neighborNodesCloudPair(const int nodeID,
                                            const int neighborNum)
{
    std::vector<size_t> neighboursID(neighborNum);
    std::vector<double> outDistSqr(neighborNum);
    kdTree_.query(nodeID, neighborNum, &neighboursID[0], &outDistSqr[0]);

    nodesCloud cloud(neighborNum);
    std::sort(neighboursID.begin(), neighboursID.end());
    for (int i = 0; i < neighborNum; i++)
    {
        cloud.id[i] = neighboursID[i];
        cloud.nodes[i] = nodes_[neighboursID[i]];
    }
    return cloud;
}

std::vector<vec3d<double>>& MeshData::nodes()
{
    return nodes_;
}

std::shared_ptr<BoundaryCondition> MeshData::nodeBC(const int nodeID) const
{
    return nodesToBC_[nodeID];
}

int MeshData::numOfNodes() const
{
    return numOfNodes_;
};
