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

    // std::map<size_t, std::string> nodeToGoupMap_;
    // std::map<std::string, boundaryCondition> groupToBoundaryConditionMap_;
    // std::map<size_t, boundaryCondition> nodeToBoundaryConditionMap_;

    const std::vector<vec3d<double>>& nodes() const;
    const vec3d<double>& node(const size_t nodeID) const;
    const nodesCloud& nodesCloudByID(const size_t nodeID) const;
    size_t numOfNodes() const;
    std::shared_ptr<BoundaryCondition> nodeBC(const size_t nodeID) const;

    // std::shared_ptr<BoundaryCondition> nodesBC(const size_t nodeID)
    // const;

   private:
    std::shared_ptr<controlData> controlData_;
    std::vector<vec3d<double>> nodes_;
    std::vector<vec3d<double>> normals_;

    std::map<std::string, std::vector<size_t>> groupToNodesMap_;
    std::map<std::string, std::shared_ptr<BoundaryCondition>> groupToBCMap_;
    std::vector<std::shared_ptr<BoundaryCondition>> nodesToBC_;
    std::vector<std::string> nodesToGroup_;
    std::vector<nodesCloud> nodesCloud_;
    size_t numOfNodes_;
    KDTreeEigenAdaptor<std::vector<vec3d<double>>, double, 3> kdTree_;

    void compactGroupToNodesMap(
        const std::map<std::string, std::vector<size_t>>&
            groupToNodesMapBeforeCompact);

    void buildBoundaryConditions();
    void buildNodeClouds();
};

#endif
