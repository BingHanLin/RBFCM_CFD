#ifndef SCALARCONDITIONPOOL_HPP
#define SCALARCONDITIONPOOL_HPP

#include "boundaryCondition.hpp"
#include "controlData.hpp"
#include "initialCondition.hpp"
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

    InitialCondition* IC() const;
    BoundaryCondition* BCByNodeID(const size_t nodeID) const;

   private:
    ControlData* controlData_;
    MeshData* meshData_;
    std::unique_ptr<InitialCondition> IC_;
    std::map<std::string, std::unique_ptr<BoundaryCondition>> groupToBCMap_;
    std::vector<std::string> nodesToBCName_;

    void buildInitialConditions();
    void buildBoundaryConditions();
    void buildNodesToConditions();

    std::vector<std::string> BCNames() const;
};

#endif
