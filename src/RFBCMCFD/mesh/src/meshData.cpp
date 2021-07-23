#include "meshData.hpp"
#include "KDTreeEigenAdaptor.hpp"
#include "constant.hpp"
#include "messages.hpp"
#include "readFromMsh.hpp"
#include "rectangle.hpp"

#include <iostream>

MeshData::MeshData(ControlData* controlData)
    : controlData_(controlData),

      nodes_(),
      normals_(),
      nodesToGroupName_(),
      groupNameToNodesMap_()
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
    buildNodeToGroupName(groupToNodesMapNotCompact);
}

void MeshData::buildNodeClouds()
{
    const auto kdTree =
        KDTreeEigenAdaptor<std::vector<vec3d<double>>, double, 3>(nodes_);

    const double neighborRadius =
        controlData_->paramsDataAt({"solverControl", "neighborRadius"});

    clouds_.resize(nodes_.size());
    for (size_t nodeID = 0; nodeID < nodes_.size(); nodeID++)
    {
        std::vector<size_t> neighboursID;
        kdTree.query(nodeID, neighborRadius, neighboursID);

        std::sort(neighboursID.begin(), neighboursID.end());

        nodesCloud cloud(nodeID, neighboursID);
        clouds_[nodeID] = cloud;
    }
}

void MeshData::buildNodeToGroupName(
    const std::map<std::string, std::vector<size_t>>& groupToNodesMapNotCompact)
{
    nodesToGroupName_.resize(nodes_.size());
    std::fill(nodesToGroupName_.begin(), nodesToGroupName_.end(),
              NOTDEFINEDGROUPNAME);

    const auto names = controlData_->groupNames();

    for (const auto& [groupName, ids] : groupToNodesMapNotCompact)
    {
        if (std::find(names.begin(), names.end(), groupName) == names.end())
            continue;

        // TODO: parallel
        for (const size_t& id : ids)
        {
            nodesToGroupName_[id] = groupName;
        }
    }

    for (size_t i = 0; i < nodesToGroupName_.size(); i++)
    {
        groupNameToNodesMap_[nodesToGroupName_[i]].push_back(i);
    }
}

const vec3d<double>& MeshData::nodeByID(const size_t nodeID) const
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

const nodesCloud& MeshData::cloudByID(const size_t nodeID) const
{
    if (nodeID < clouds_.size())
    {
        return clouds_[nodeID];
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

const std::string& MeshData::groupNameByID(const size_t nodeID) const
{
    if (nodeID < nodesToGroupName_.size())
    {
        return nodesToGroupName_[nodeID];
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
