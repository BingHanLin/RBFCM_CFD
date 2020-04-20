#ifndef BOUNDARYCONDITIONS_HPP
#define BOUNDARYCONDITIONS_HPP

#include "MQBasis.hpp"
#include "enumMap.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>


class MeshData;

class BoundaryCondition
{
   public:
    BoundaryCondition(){};

    virtual boundaryConditionType type() = 0;
    //  virtual void setRHSValue(){};
    virtual void fillCoeffMatrix(const int nodeID,
                                 const std::vector<int>& nodesCloudID,
                                 const std::vector<vec3d<double>>& nodesCloud,
                                 std::shared_ptr<MQBasis> RBFBasis,
                                 Eigen::SparseMatrix<double>& spMatrix) = 0;

   protected:
    MeshData* mesh_;
};

#endif
