#ifndef PUCONDITIONPOOL_HPP
#define PUCONDITIONPOOL_HPP

#include "BoundaryCondition.hpp"
#include "ControlData.hpp"
#include "InitialCondition.hpp"
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
class PUConditionPool
{
   public:
    explicit PUConditionPool(ControlData* controlData, MeshData* meshData);

    InitialCondition* UICByNodeID(const size_t nodeID) const;
    BoundaryCondition* UBCByNodeID(const size_t nodeID) const;

    InitialCondition* PICByNodeID(const size_t nodeID) const;
    BoundaryCondition* PBCByNodeID(const size_t nodeID) const;

   private:
    ControlData* controlData_;
    MeshData* meshData_;
    std::map<std::string, std::unique_ptr<InitialCondition>> groupToUICMap_;
    std::map<std::string, std::unique_ptr<BoundaryCondition>> groupToUBCMap_;

    std::map<std::string, std::unique_ptr<InitialCondition>> groupToPICMap_;
    std::map<std::string, std::unique_ptr<BoundaryCondition>> groupToPBCMap_;

    void buildInitialConditions();
    void buildBoundaryConditions();
};

#endif
