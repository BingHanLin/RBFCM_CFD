#ifndef MESHDATA_HPP
#define MESHDATA_HPP

#include "dataStructure.hpp"
#include "enumMap.hpp"
#include "json.h"
#include "readFromMsh.hpp"
#include <vector>


/*************************************************************************


*************************************************************************/

class MeshData
{
   public:
    MeshData(nlohmann::json& geometryControlParams);
    ~MeshData(){};

    nlohmann::json& geometryControlParams_;
    int numInnNodes_, numBouNodes_, numAllNodes_;
    std::map<std::string, std::vector<int>> groupToNodesMap_;
    std::map<std::string, boundaryCondition> groupToBoundaryConditionMap_;
    std::map<int, std::string> nodeToGoupMap_;
    std::map<int, boundaryCondition> nodeToBoundaryConditionMap_;

    std::vector<vec3d<double>> nodes_;
    std::vector<vec3d<double>> normals_;

   private:
    std::string meshFilePath_;

    void buildMesh();
    void readMeshFromFile();
};

#endif
