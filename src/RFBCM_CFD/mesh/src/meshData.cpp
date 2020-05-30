#include "meshData.hpp"
#include "constantValueBC.hpp"
#include "readFromMsh.hpp"
#include "rectangle.hpp"

#include <iostream>

MeshData::MeshData(std::shared_ptr<controlData> inControlData)
    : controlData_(inControlData),
      numOfNodes_(),
      groupToNodesMap_(),
      nodes_(),
      normals_(),
      nodesToGroup_()
{
    std::cout << "MeshData" << std::endl;

    meshType meshType = controlData_->paramsDataAt({"geometryControl", "type"});

    std::map<std::string, std::vector<int>> groupToNodesMapBeforeCompact;

    if (meshType == meshType::DEFAULT)
    {
        std::cout << "read nodes from msh file" << std::endl;

        const std::string meshFileName =
            controlData_->paramsDataAt({"geometryControl", "fileName"})
                .get<std::string>();

        const std::string absPath =
            controlData_->workingDir().concat("/" + meshFileName).string();

        bool isReadSuccess = readFromMsh(absPath, nodes_, normals_,
                                         groupToNodesMapBeforeCompact);
    }
    else if (meshType == meshType::RECTNAGLE)
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
    compactGroupToNodesMap(groupToNodesMapBeforeCompact);
    buildBoundaryConditions();
}

void MeshData::compactGroupToNodesMap(
    const std::map<std::string, std::vector<int>>& groupToNodesMapBeforeCompact)
{
    nodesToGroup_.resize(nodes_.size());
    groupToNodesMap_.clear();

    const auto& constantValueBCData = controlData_->paramsDataAt(
        {"physicsControl", "boundaryConditions", "constantValue"});

    for (const auto& oneBCData : constantValueBCData)
    {
        const std::string groupName = oneBCData.at("groupName");
        groupToNodesMap_.insert({groupName, {}});

        // TODO: parallel
        for (const int& nodeIndex : groupToNodesMapBeforeCompact.at(groupName))
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
        {"physicsControl", "boundaryConditions", "constantValue"});

    for (auto& oneBCData : constantValueBCData)
    {
        const std::string groupName = oneBCData.at("groupName");

        auto constantValueBC =
            std::make_shared<ConstantValueBC>(oneBCData.at("value"));
        groupToBCMap_.insert({groupName, constantValueBC});

        for (int& nodeID : groupToNodesMap_.at(groupName))
        {
            nodesToBC_[nodeID] = constantValueBC;
        }
    }
}

nodesCloud MeshData::neighborNodesCloud(const int nodeID, const int neighborNum)
{
    std::vector<size_t> neighboursID(neighborNum);
    std::vector<double> outDistSqr(neighborNum);
    kdTree_.query(nodeID, neighborNum, &neighboursID[0], &outDistSqr[0]);

    nodesCloud cloud(neighborNum);
    for (int i = 0; i < neighborNum; i++)
    {
        cloud.id[i] = neighboursID[i];
        cloud.nodes[i] = nodes_[neighboursID[i]];
    }

    // std::sort(neighboursID.begin(), neighboursID.end());
    // for (int i = 0; i < neighborNum; i++)
    // {
    //     cloud.id[i] = neighboursID[i];
    //     cloud.nodes[i] = nodes_[neighboursID[i]];
    // }
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
