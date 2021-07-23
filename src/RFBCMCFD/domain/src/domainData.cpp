#include "domainData.hpp"
#include "boundaryCondition.hpp"
#include "constant.hpp"
#include "constantValueBC.hpp"
#include "meshData.hpp"
#include "messages.hpp"
#include "neumannBC.hpp"

#include "constantValueIC.hpp"
#include "initialCondition.hpp"

#include <iostream>

DomainData::DomainData(ControlData* controlData, MeshData* meshData)
    : controlData_(controlData), meshData_(meshData)
{
    buildInitialConditions();
    buildBoundaryConditions();
}

DomainData::~DomainData(){};

void DomainData::buildInitialConditions()
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

void DomainData::buildBoundaryConditions()
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

BoundaryCondition* DomainData::BCByID(const size_t nodeID) const
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

InitialCondition* DomainData::ICByID(const size_t nodeID) const
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

MeshData* DomainData::meshData() const
{
    return meshData_;
}
