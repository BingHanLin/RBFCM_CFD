#ifndef MESHDATA_HPP
#define MESHDATA_HPP

#include "BoundaryCondition.hpp"
#include "ControlData.hpp"
#include "dataStructure.hpp"
#include "enumMap.hpp"
#include "json.h"
#include <vector>

/*************************************************************************


*************************************************************************/

class MeshData
{
   public:
    MeshData(ControlData* controlData);
    ~MeshData(){};

    const std::vector<vec3d<double>>& nodes() const;
    size_t numOfNodes() const;

    const vec3d<double>& nodeByID(const size_t nodeID) const;
    const nodesCloud& cloudByID(const size_t nodeID) const;
    const vec3d<double>& normalByID(const size_t nodeID) const;

    const std::map<std::string, std::vector<size_t>>& groupNameToNodesMap()
        const;

   private:
    ControlData* controlData_;
    std::vector<vec3d<double>> nodes_;
    std::vector<vec3d<double>> normals_;
    std::vector<nodesCloud> clouds_;
    std::map<std::string, std::vector<size_t>>
        groupNameToNodesMap_;  // may not compact

    void buildNodeClouds();
};

#endif
