#ifndef MESHDATA_HPP
#define MESHDATA_HPP
#include "boundaryCondition.hpp"

#include "dataStructure.hpp"
#include "enumMap.hpp"
#include "json.h"
#include <vector>

/*************************************************************************


*************************************************************************/

class MeshData
{
   public:
    MeshData(nlohmann::json& controlParams);
    ~MeshData(){};

    nlohmann::json& geometryControlParams_;
    nlohmann::json& physicsControlParams_;
    // std::map<int, std::string> nodeToGoupMap_;
    // std::map<std::string, boundaryCondition> groupToBoundaryConditionMap_;
    // std::map<int, boundaryCondition> nodeToBoundaryConditionMap_;

    int numOfNodes() const;

   private:
    std::string meshFilePath_;
    std::vector<vec3d<double>> nodes_;
    std::vector<vec3d<double>> normals_;
    std::map<std::string, std::vector<int>> groupToNodesMap_;
    std::map<std::string, std::shared_ptr<BoundaryCondition>> groupToBCMap_;
    std::vector<std::shared_ptr<BoundaryCondition>> nodesBC_;

    int numOfNodes_;

    // void compactGroupToNodesMap(
    //     std::map<std::string, std::vector<int>>&
    //     groupToNodesMapBeforeCompact);
    void buildBoundaryConditions();
    void matchBoundaryConditions(
        std::map<std::string, std::vector<int>> groupToNodesMapBeforeCompact);
};

#endif
