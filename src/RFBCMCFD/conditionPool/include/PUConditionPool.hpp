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

    InitialCondition* UIC() const;
    BoundaryCondition* UBCByNodeID(const size_t nodeID) const;

    InitialCondition* PIC() const;
    BoundaryCondition* PBCByNodeID(const size_t nodeID) const;

   private:
    ControlData* controlData_;
    MeshData* meshData_;
    std::unique_ptr<InitialCondition> UIC_;
    std::map<std::string, std::unique_ptr<BoundaryCondition>> groupToUBCMap_;
    std::vector<std::string> nodesToUBCName_;

    std::unique_ptr<InitialCondition> PIC_;
    std::map<std::string, std::unique_ptr<BoundaryCondition>> groupToPBCMap_;
    std::vector<std::string> nodesToPBCName_;

    void buildInitialConditions();
    void buildBoundaryConditions();
    void buildNodesToConditions();

    std::vector<std::string> BCNames() const;
};

#endif
