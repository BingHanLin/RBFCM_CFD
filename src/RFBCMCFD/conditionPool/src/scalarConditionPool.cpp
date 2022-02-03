#include "ScalarConditionPool.hpp"
#include "MeshData.hpp"
#include "constant.hpp"
#include "messages.hpp"

#include "ConstantValueBC.hpp"
#include "NeumannBC.hpp"

#include "ConstantValueIC.hpp"
#include "InitialCondition.hpp"

#include <iostream>

ScalarConditionPool::ScalarConditionPool(ControlData* controlData,
                                         MeshData* meshData)
    : controlData_(controlData), meshData_(meshData)
{
    buildInitialConditions();
    buildBoundaryConditions();
}

void ScalarConditionPool::buildInitialConditions()
{
    const auto& constValueICData =
        controlData_->paramsDataAt({"initialConditions", "constantValue"});

    IC_ = std::make_unique<ConstantValueIC>(constValueICData.at("value"));
}
void ScalarConditionPool::buildBoundaryConditions()
{
    const auto& NeumannBCData =
        controlData_->paramsDataAt({"boundaryConditions", "neumann"});

    for (auto& oneBCData : NeumannBCData)
    {
        const std::string groupName = oneBCData.at("groupName");

        auto bc = std::make_unique<NeumannBC>(oneBCData.at("value"), meshData_);

        groupToBCMap_.insert({groupName, std::move(bc)});
    }

    const auto& constantValueBCData =
        controlData_->paramsDataAt({"boundaryConditions", "constantValue"});

    for (auto& oneBCData : constantValueBCData)
    {
        const std::string groupName = oneBCData.at("groupName");

        auto bc =
            std::make_unique<ConstantValueBC>(oneBCData.at("value"), meshData_);
        groupToBCMap_.insert({groupName, std::move(bc)});
    }
}

void ScalarConditionPool::buildNodesToConditions()
{
    const auto numOfNodes = meshData_->numOfNodes();
    const auto names = this->BCNames();

    nodesToBCName_.resize(numOfNodes);
    std::fill(nodesToBCName_.begin(), nodesToBCName_.end(),
              NOTDEFINED_GROUPNAME);

    const auto groupNameToNodesMap = meshData_->groupNameToNodesMap();

    for (const auto& oneBCName : this->BCNames())
    {
        if (groupNameToNodesMap.find(oneBCName) == groupNameToNodesMap.end())
            continue;

        for (const size_t& id : groupNameToNodesMap.at(oneBCName))
            nodesToBCName_[id] = oneBCName;
    }
}

InitialCondition* ScalarConditionPool::IC() const
{
    return IC_.get();
}

BoundaryCondition* ScalarConditionPool::BCByNodeID(const size_t nodeID) const
{
    if (groupToBCMap_.find(nodesToBCName_[nodeID]) != groupToBCMap_.end())
    {
        return groupToBCMap_.at(nodesToBCName_[nodeID]).get();
    }
    else
    {
        return nullptr;
    }
}

std::vector<std::string> ScalarConditionPool::BCNames() const
{
    std::vector<std::string> result;

    const auto& NeumannBCData =
        controlData_->paramsDataAt({"boundaryConditions", "neumann"});

    for (auto& oneBCData : NeumannBCData)
    {
        const std::string groupName = oneBCData.at("groupName");
        result.push_back(groupName);
    }

    const auto& constantValueBCData =
        controlData_->paramsDataAt({"boundaryConditions", "constantValue"});

    for (auto& oneBCData : constantValueBCData)
    {
        const std::string groupName = oneBCData.at("groupName");
        result.push_back(groupName);
    }

    return result;
}
