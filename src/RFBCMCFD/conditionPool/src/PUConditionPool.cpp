#include "PUConditionPool.hpp"
#include "MeshData.hpp"
#include "constant.hpp"
#include "messages.hpp"

#include "ConstantValueBC.hpp"
#include "NeumannBC.hpp"

#include "ConstantValueIC.hpp"
#include "InitialCondition.hpp"

#include <iostream>

PUConditionPool::PUConditionPool(ControlData* controlData, MeshData* meshData)
    : controlData_(controlData), meshData_(meshData)
{
    buildInitialConditions();
    buildBoundaryConditions();
}

void PUConditionPool::buildInitialConditions()
{
    const auto& constValueUICData =
        controlData_->paramsDataAt({"initialConditions", "U", "constantValue"});

    UIC_ = std::make_unique<ConstantValueIC>(constValueUICData.at("value"));
}

void PUConditionPool::buildBoundaryConditions()
{
    const auto& neumannUBCData =
        controlData_->paramsDataAt({"boundaryConditions", "U", "neumann"});

    for (auto& oneBCData : neumannUBCData)
    {
        const std::string groupName = oneBCData.at("groupName");

        auto bc = std::make_unique<neumannBC>(oneBCData.at("value"), meshData_);

        groupToUBCMap_.insert({groupName, std::move(bc)});
    }

    const auto& constantValueUBCData = controlData_->paramsDataAt(
        {"boundaryConditions", "U", "constantValue"});

    for (auto& oneBCData : constantValueUBCData)
    {
        const std::string groupName = oneBCData.at("groupName");

        auto bc =
            std::make_unique<ConstantValueBC>(oneBCData.at("value"), meshData_);
        groupToUBCMap_.insert({groupName, std::move(bc)});
    }
}

InitialCondition* PUConditionPool::UIC() const
{
    return UIC_.get();
}

BoundaryCondition* PUConditionPool::UBCByNodeID(const size_t nodeID) const
{
    if (groupToUBCMap_.find(meshData_->groupNameByID(nodeID)) !=
        groupToUBCMap_.end())
    {
        return groupToUBCMap_.at(meshData_->groupNameByID(nodeID)).get();
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
    if (groupToPBCMap_.find(meshData_->groupNameByID(nodeID)) !=
        groupToPBCMap_.end())
    {
        return groupToPBCMap_.at(meshData_->groupNameByID(nodeID)).get();
    }
    else
    {
        return nullptr;
    }
}
