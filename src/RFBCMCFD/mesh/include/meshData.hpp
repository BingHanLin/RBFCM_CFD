#ifndef MESHDATA_HPP
#define MESHDATA_HPP

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
    MeshData(ControlData* controlData);
    ~MeshData(){};

    const std::vector<vec3d<double>>& nodes() const;
    size_t numOfNodes() const;

    const vec3d<double>& nodeByID(const size_t nodeID) const;
    const nodesCloud& cloudByID(const size_t nodeID) const;
    const vec3d<double>& normalByID(const size_t nodeID) const;
    const std::string& groupNameByID(const size_t nodeID) const;

   private:
    ControlData* controlData_;
    std::vector<vec3d<double>> nodes_;
    std::vector<vec3d<double>> normals_;
    std::vector<nodesCloud> clouds_;
    std::vector<std::string> nodesToGroupName_;

    void buildNodeToGroupName(const std::map<std::string, std::vector<size_t>>&
                                  groupToNodesMapNotCompact);
    void buildNodeClouds();
};

#endif
