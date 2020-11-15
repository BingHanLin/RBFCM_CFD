#include "meshData.hpp"
#include "constant.hpp"
#include "constantValueBC.hpp"
#include "messages.hpp"
#include "neumannBC.hpp"
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

void MeshData::compactGroupToNodesMap(
    const std::map<std::string, std::vector<size_t>>& groupToNodesMapNotCompact)
{
    groupToNodesMap_.clear();

    nodesToGroup_.resize(nodes_.size());
    std::fill(nodesToGroup_.begin(), nodesToGroup_.end(), NOTDEFINEDGROUPNAME);

    const auto names = controlData_->BCGroupNames();

    for (const auto& [groupName, ids] : groupToNodesMapNotCompact)
    {
        if (std::find(names.begin(), names.end(), groupName) == names.end())
            continue;

        groupToNodesMap_.insert({groupName, {}});

        // TODO: parallel
        for (const size_t& id : ids)
        {
            nodesToGroup_[id] = groupName;
        }
    }

    for (size_t id = 0; id < nodesToGroup_.size(); id++)
    {
        groupToNodesMap_[nodesToGroup_[id]].push_back(id);
    }
}

void MeshData::buildBoundaryConditions()
{
    nodesToBC_.resize(nodes_.size());

    const auto& neumannBCData = controlData_->paramsDataAt(
        {"physicsControl", "boundaryConditions", "neumann"});

    for (auto& oneBCData : neumannBCData)
    {
        const std::string groupName = oneBCData.at("groupName");

        if (groupToNodesMap_.find(groupName) == groupToNodesMap_.end())
        {
            ASSERT("buildBoundaryConditions: "
                   << groupName << "not found in groupToNodesMap_");
        }

        auto bc = std::make_shared<neumannBC>(oneBCData.at("value"), this);
        groupToBCMap_.insert({groupName, bc});

        for (size_t& nodeID : groupToNodesMap_.at(groupName))
        {
            nodesToBC_[nodeID] = bc;
        }
    }

    const auto& constantValueBCData = controlData_->paramsDataAt(
        {"physicsControl", "boundaryConditions", "constantValue"});

    for (auto& oneBCData : constantValueBCData)
    {
        const std::string groupName = oneBCData.at("groupName");

        if (groupToNodesMap_.find(groupName) == groupToNodesMap_.end())
        {
            ASSERT("buildBoundaryConditions: "
                   << groupName << "not found in groupToNodesMap_");
        }

        auto bc =
            std::make_shared<ConstantValueBC>(oneBCData.at("value"), this);
        groupToBCMap_.insert({groupName, bc});

        for (size_t& nodeID : groupToNodesMap_.at(groupName))
        {
            nodesToBC_[nodeID] = bc;
        }
    }
}

const std::vector<size_t>& MeshData::nodesIDByGroupName(
    const std::string groupName) const
{
    return groupToNodesMap_.at(groupName);
}

const nodesCloud& MeshData::nodesCloudByID(const size_t nodeID) const
{
    if (nodeID < nodesCloud_.size())
    {
        return nodesCloud_[nodeID];
    }
    else
    {
        ASSERT("node ID out of bound, ID: " << nodeID);
    }
}

const vec3d<double>& MeshData::normalByID(const size_t nodeID) const
{
    if (nodeID < normals_.size())
    {
        return normals_[nodeID];
    }
    else
    {
        ASSERT("node ID out of bound, ID: " << nodeID);
    }
}

// const vec3d<double>& MeshData::node(const size_t nodeID) const
// {
//     if (nodeID < nodes_.size())
//     {
//         return nodes_[nodeID];
//     }
//     else
//     {
//         ASSERT("node ID out of bound, ID: " << nodeID);
//     }
// }

std::shared_ptr<BoundaryCondition> MeshData::nodeBCByID(
    const size_t nodeID) const
{
    if (nodeID < nodesToBC_.size())
    {
        return nodesToBC_[nodeID];
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

size_t MeshData::numOfNodes() const
{
    return nodes_.size();
};
