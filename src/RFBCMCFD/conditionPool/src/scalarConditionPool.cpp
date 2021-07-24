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

    for (auto& oneICData : constValueICData)
    {
        const std::string groupName = oneICData.at("groupName");

        auto ic = std::make_unique<ConstantValueIC>(oneICData.at("value"));

        groupToICMap_.insert({groupName, std::move(ic)});
    }
}
void ScalarConditionPool::buildBoundaryConditions()
{
    const auto& neumannBCData =
        controlData_->paramsDataAt({"boundaryConditions", "neumann"});

    for (auto& oneBCData : neumannBCData)
    {
        const std::string groupName = oneBCData.at("groupName");

        auto bc = std::make_unique<neumannBC>(oneBCData.at("value"), meshData_);

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

InitialCondition* ScalarConditionPool::ICByNodeID(const size_t nodeID) const
{
    if (groupToICMap_.find(meshData_->groupNameByID(nodeID)) !=
        groupToICMap_.end())
    {
        return groupToICMap_.at(meshData_->groupNameByID(nodeID)).get();
    }
    else
    {
        return nullptr;
    }
}

BoundaryCondition* ScalarConditionPool::BCByNodeID(const size_t nodeID) const
{
    if (groupToBCMap_.find(meshData_->groupNameByID(nodeID)) !=
        groupToBCMap_.end())
    {
        return groupToBCMap_.at(meshData_->groupNameByID(nodeID)).get();
    }
    else
    {
        return nullptr;
    }
}
