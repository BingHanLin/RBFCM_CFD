#include "scalarBCPool.hpp"
#include "constant.hpp"
#include "constantValueBC.hpp"
#include "meshData.hpp"
#include "messages.hpp"
#include "neumannBC.hpp"

#include <iostream>

ScalarBCPool::ScalarBCPool(ControlData* controlData, MeshData* meshData)
    : controlData_(controlData), meshData_(meshData)
{
    buildBoundaryConditions();
}

void ScalarBCPool::buildBoundaryConditions()
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

BoundaryCondition* ScalarBCPool::BCByNodeID(const size_t nodeID) const
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
