#include "meshData.hpp"
#include "constantValueBC.hpp"
#include "messages.hpp"
#include "readFromMsh.hpp"
#include "rectangle.hpp"

#include <iostream>

MeshData::MeshData(std::shared_ptr<controlData> inControlData)
    : controlData_(inControlData),
      groupToNodesMap_(),
      nodes_(),
      normals_(),
      nodesToGroup_()
{
    std::cout << "MeshData" << std::endl;

    meshType meshType = controlData_->paramsDataAt({"geometryControl", "type"});

    std::map<std::string, std::vector<size_t>> groupToNodesMapNotCompact;

    if (meshType == meshType::DEFAULT)
    {
        std::cout << "read nodes from msh file" << std::endl;

        const std::string meshFileName =
            controlData_->paramsDataAt({"geometryControl", "fileName"})
                .get<std::string>();

        const std::string absPath =
            controlData_->workingDir().concat("/" + meshFileName).string();

        bool isReadSuccess =
            readFromMsh(absPath, nodes_, normals_, groupToNodesMapNotCompact);
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

    buildNodeClouds();
    compactGroupToNodesMap(groupToNodesMapNotCompact);
    buildBoundaryConditions();
}

void MeshData::compactGroupToNodesMap(
    const std::map<std::string, std::vector<size_t>>& groupToNodesMapNotCompact)
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
        for (const size_t& nodeIndex : groupToNodesMapNotCompact.at(groupName))
        {
            nodesToGroup_[nodeIndex] = groupName;
        }
    }

    for (size_t nodeID = 0; nodeID < nodes_.size(); nodeID++)
    {
        groupToNodesMap_[nodesToGroup_[nodeID]].push_back(nodeID);
    }
}

void MeshData::buildNodeClouds()
{
    kdTree_ = KDTreeEigenAdaptor<std::vector<vec3d<double>>, double, 3>(nodes_);

    size_t neighborNum =
        controlData_->paramsDataAt({"solverControl", "neighborNumber"});

    nodesCloud_.resize(nodes_.size());
    for (size_t nodeID = 0; nodeID < nodes_.size(); nodeID++)
    {
        std::vector<size_t> neighboursID(neighborNum);
        std::vector<double> outDistSqr(neighborNum);

        kdTree_.query(nodeID, neighborNum, &neighboursID[0], &outDistSqr[0]);

        std::sort(neighboursID.begin(), neighboursID.end());

        nodesCloud cloud(nodeID, neighboursID);
        nodesCloud_[nodeID] = cloud;
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
            std::make_shared<ConstantValueBC>(oneBCData.at("value"), this);
        groupToBCMap_.insert({groupName, constantValueBC});

        for (size_t& nodeID : groupToNodesMap_.at(groupName))
        {
            nodesToBC_[nodeID] = constantValueBC;
        }
    }
}

const nodesCloud& MeshData::nodesCloudByID(const size_t nodeID) const
{
    return nodesCloud_[nodeID];
}

const vec3d<double>& MeshData::node(const size_t nodeID) const
{
    if (nodeID < nodes_.size())
    {
        return nodes_[nodeID];
    }
    else
    {
        ASSERT("node ID out of bound, ID: " << nodeID);
    }
}

const std::vector<vec3d<double>>& MeshData::nodes() const
{
    return nodes_;
}

std::shared_ptr<BoundaryCondition> MeshData::nodeBC(const size_t nodeID) const
{
    return nodesToBC_[nodeID];
}

size_t MeshData::numOfNodes() const
{
    return nodes_.size();
};
