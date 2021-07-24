#ifndef SCALARBCPOOL_HPP
#define SCALARBCPOOL_HPP

#include "boundaryCondition.hpp"
#include "controlData.hpp"
#include "dataStructure.hpp"
#include "enumMap.hpp"
#include "json.h"
#include <memory>
#include <vector>

/*************************************************************************


*************************************************************************/
class BoundaryCondition;
class MeshData;
class ScalarBCPool
{
   public:
    explicit ScalarBCPool(ControlData* controlData, MeshData* meshData);

    BoundaryCondition* BCByNodeID(const size_t nodeID) const;

   private:
    ControlData* controlData_;
    MeshData* meshData_;
    std::map<std::string, std::unique_ptr<BoundaryCondition>> groupToBCMap_;

    void buildBoundaryConditions();
};

#endif
