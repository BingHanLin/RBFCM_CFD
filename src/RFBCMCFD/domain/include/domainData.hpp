#ifndef DOMAINDATA_HPP
#define DOMAINDATA_HPP

#include "controlData.hpp"
#include "dataStructure.hpp"
#include "enumMap.hpp"
#include "json.h"
#include <memory>
#include <vector>

/*************************************************************************


*************************************************************************/
class BoundaryCondition;
class InitialCondition;
class MeshData;
class DomainData
{
   public:
    DomainData(ControlData* controlData, MeshData* meshData);
    ~DomainData();

    MeshData* meshData() const;
    BoundaryCondition* BCByID(const size_t nodeID) const;
    InitialCondition* ICByID(const size_t nodeID) const;

   private:
    ControlData* controlData_;
    MeshData* meshData_;
    std::map<std::string, std::unique_ptr<BoundaryCondition>> groupToBCMap_;
    std::map<std::string, std::unique_ptr<InitialCondition>> groupToICMap_;

    void buildInitialConditions();
    void buildBoundaryConditions();
};

#endif
