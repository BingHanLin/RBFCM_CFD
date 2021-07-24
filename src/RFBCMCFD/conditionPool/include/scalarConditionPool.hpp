#ifndef SCALARCONDITIONPOOL_HPP
#define SCALARCONDITIONPOOL_HPP

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
class ScalarConditionPool
{
   public:
    explicit ScalarConditionPool(ControlData* controlData, MeshData* meshData);

    InitialCondition* ICByNodeID(const size_t nodeID) const;
    BoundaryCondition* BCByNodeID(const size_t nodeID) const;

   private:
    ControlData* controlData_;
    MeshData* meshData_;
    std::map<std::string, std::unique_ptr<InitialCondition>> groupToICMap_;
    std::map<std::string, std::unique_ptr<BoundaryCondition>> groupToBCMap_;

    void buildInitialConditions();
    void buildBoundaryConditions();
};

#endif
