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

struct nodesCloud
{
    nodesCloud(const int size)
    {
        id.resize(size);
        nodes.resize(size);
    }
    std::vector<int> id;
    std::vector<vec3d<double>> nodes;
};

class MeshData
{
   public:
    MeshData();
    ~MeshData(){};

    // std::map<int, std::string> nodeToGoupMap_;
    // std::map<std::string, boundaryCondition> groupToBoundaryConditionMap_;
    // std::map<int, boundaryCondition> nodeToBoundaryConditionMap_;

    std::vector<vec3d<double>>& nodes();
    int numOfNodes() const;
    std::shared_ptr<BoundaryCondition> nodeBC(const int nodeID) const;

    nodesCloud neighborNodesCloudPair(const int nodeID, const int neighborNum);

    // std::shared_ptr<BoundaryCondition> nodesBC(const int nodeID) const;

   private:
    controlData* controlData_;
    std::vector<vec3d<double>> nodes_;
    std::vector<vec3d<double>> normals_;

    std::map<std::string, std::vector<int>> groupToNodesMap_;
    std::map<std::string, std::shared_ptr<BoundaryCondition>> groupToBCMap_;
    std::vector<std::shared_ptr<BoundaryCondition>> nodesToBC_;
    std::vector<std::string> nodesToGroup_;

    int numOfNodes_;
    KDTreeEigenAdaptor<std::vector<vec3d<double>>, double, 3> kdTree_;

    void buildBoundaryConditions();
    void compactGroupToNodesMap(
        std::map<std::string, std::vector<int>> groupToNodesMapBeforeCompact);
};

#endif
