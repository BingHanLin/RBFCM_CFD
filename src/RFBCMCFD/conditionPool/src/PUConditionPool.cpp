#include "PUConditionPool.hpp"
#include "MeshData.hpp"
#include "constant.hpp"
#include "messages.hpp"

#include "ConstantVecValueBC.hpp"

#include "ConstantVecValueIC.hpp"

#include <iostream>

PUConditionPool::PUConditionPool(ControlData* controlData, MeshData* meshData)
    : controlData_(controlData), meshData_(meshData)
{
    buildInitialConditions();
    buildBoundaryConditions();
    buildNodesToConditions();
}

void PUConditionPool::buildInitialConditions()
{
    const auto& constValueUICData =
        controlData_->paramsDataAt({"initialConditions", "U", "constantValue"});

    const auto vel = constValueUICData.at("value").get<std::array<double, 3>>();

    UIC_ = std::make_unique<ConstantVecValueIC>(vel, meshData_);
}

void PUConditionPool::buildBoundaryConditions()
{
    const auto& constantValueUBCData = controlData_->paramsDataAt(
        {"boundaryConditions", "U", "constantValue"});

    for (auto& oneBCData : constantValueUBCData)
    {
        const std::string groupName = oneBCData.at("groupName");

        const auto vel = oneBCData.at("value").get<std::array<double, 3>>();

        auto bc = std::make_unique<ConstantVecValueBC>(vel, meshData_);
        groupToUBCMap_.insert({groupName, std::move(bc)});
    }
}

void PUConditionPool::buildNodesToConditions()
{
    const auto numOfNodes = meshData_->numOfNodes();
    const auto names = this->BCNames();

    nodesToUBCName_.resize(numOfNodes);
    std::fill(nodesToUBCName_.begin(), nodesToUBCName_.end(),
              NOTDEFINED_GROUPNAME);

    const auto groupNameToNodesMap = meshData_->groupNameToNodesMap();

    for (const auto& oneBCName : this->BCNames())
    {
        if (groupNameToNodesMap.find(oneBCName) == groupNameToNodesMap.end())
            continue;

        for (const size_t& id : groupNameToNodesMap.at(oneBCName))
            nodesToUBCName_[id] = oneBCName;
    }
}

InitialCondition* PUConditionPool::UIC() const
{
    return UIC_.get();
}

BoundaryCondition* PUConditionPool::UBCByNodeID(const size_t nodeID) const
{
    if (groupToUBCMap_.find(nodesToUBCName_[nodeID]) != groupToUBCMap_.end())
    {
        return groupToUBCMap_.at(nodesToUBCName_[nodeID]).get();
    }
    else
    {
        return nullptr;
    }
}

InitialCondition* PUConditionPool::PIC() const
{
    return PIC_.get();
}

BoundaryCondition* PUConditionPool::PBCByNodeID(const size_t nodeID) const
{
    if (groupToPBCMap_.find(nodesToPBCName_[nodeID]) != groupToPBCMap_.end())
    {
        return groupToPBCMap_.at(nodesToPBCName_[nodeID]).get();
    }
    else
    {
        return nullptr;
    }
}

std::vector<std::string> PUConditionPool::BCNames() const
{
    std::vector<std::string> result;

    const auto& constantValueUBCData = controlData_->paramsDataAt(
        {"boundaryConditions", "U", "constantValue"});

    for (auto& oneBCData : constantValueUBCData)
    {
        const std::string groupName = oneBCData.at("groupName");
        result.push_back(groupName);
    }

    return result;
}
