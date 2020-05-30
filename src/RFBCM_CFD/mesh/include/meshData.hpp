#ifndef MESHDATA_HPP
#define MESHDATA_HPP

#include "KDTreeEigenAdaptor.hpp"
#include "boundaryCondition.hpp"
#include "controlData.hpp"
#include "dataStructure.hpp"
#include "enumMap.hpp"
#include "json.h"
#include <vector>

/*************************************************************************


*************************************************************************/

class MeshData
{
   public:
    MeshData(std::shared_ptr<controlData> inControlData);
    ~MeshData(){};

    // std::map<int, std::string> nodeToGoupMap_;
    // std::map<std::string, boundaryCondition> groupToBoundaryConditionMap_;
    // std::map<int, boundaryCondition> nodeToBoundaryConditionMap_;

    std::vector<vec3d<double>>& nodes();
    int numOfNodes() const;
    std::shared_ptr<BoundaryCondition> nodeBC(const int nodeID) const;

    nodesCloud neighborNodesCloud(const int nodeID, const int neighborNum);

    // std::shared_ptr<BoundaryCondition> nodesBC(const int nodeID) const;

   private:
    std::shared_ptr<controlData> controlData_;
    std::vector<vec3d<double>> nodes_;
    std::vector<vec3d<double>> normals_;

    std::map<std::string, std::vector<int>> groupToNodesMap_;
    std::map<std::string, std::shared_ptr<BoundaryCondition>> groupToBCMap_;
    std::vector<std::shared_ptr<BoundaryCondition>> nodesToBC_;
    std::vector<std::string> nodesToGroup_;

    int numOfNodes_;
    KDTreeEigenAdaptor<std::vector<vec3d<double>>, double, 3> kdTree_;

    void buildBoundaryConditions();
    void compactGroupToNodesMap(const std::map<std::string, std::vector<int>>&
                                    groupToNodesMapBeforeCompact);
};

#endif
