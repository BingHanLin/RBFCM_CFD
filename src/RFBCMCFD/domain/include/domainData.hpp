#ifndef DOMAINDATA_HPP
#define DOMAINDATA_HPP

#include "KDTreeEigenAdaptor.hpp"
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
    DomainData(std::shared_ptr<controlData> inControlData);
    ~DomainData();

    MeshData* meshData() const;
    BoundaryCondition* BCByID(const size_t nodeID) const;
    InitialCondition* ICByID(const size_t nodeID) const;

   private:
    std::shared_ptr<controlData> controlData_;
    std::unique_ptr<MeshData> meshData_;
    std::map<std::string, std::unique_ptr<BoundaryCondition>> groupToBCMap_;
    std::map<std::string, std::unique_ptr<InitialCondition>> groupToICMap_;

    void buildInitialConditions();
    void buildBoundaryConditions();
};

#endif
