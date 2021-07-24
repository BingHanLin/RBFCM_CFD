#include "scalarICPool.hpp"
#include "constant.hpp"
#include "constantValueIC.hpp"
#include "initialCondition.hpp"
#include "meshData.hpp"
#include "messages.hpp"

#include <iostream>

ScalarICPool::ScalarICPool(ControlData* controlData, MeshData* meshData)
    : controlData_(controlData), meshData_(meshData)
{
    buildInitialConditions();
}

void ScalarICPool::buildInitialConditions()
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

InitialCondition* ScalarICPool::ICByNodeID(const size_t nodeID) const
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
