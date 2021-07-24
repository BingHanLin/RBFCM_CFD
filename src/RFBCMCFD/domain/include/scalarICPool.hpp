#ifndef SCALARICPOOL_HPP
#define SCALARICPOOL_HPP

#include "controlData.hpp"
#include "dataStructure.hpp"
#include "initialCondition.hpp"

#include "enumMap.hpp"
#include "json.h"
#include <memory>
#include <vector>

/*************************************************************************


*************************************************************************/
class InitialCondition;
class MeshData;
class ScalarICPool
{
   public:
    explicit ScalarICPool(ControlData* controlData, MeshData* meshData);

    InitialCondition* ICByNodeID(const size_t nodeID) const;

   private:
    ControlData* controlData_;
    MeshData* meshData_;
    std::map<std::string, std::unique_ptr<InitialCondition>> groupToICMap_;

    void buildInitialConditions();
};

#endif
