#ifndef BOUNDARYCONDITION_HPP
#define BOUNDARYCONDITION_HPP

#include "dataStructure.hpp"

#include <memory>

#include <Eigen/Dense>
#include <Eigen/Sparse>

class MeshData;
class MQBasis;

class BoundaryCondition
{
   public:
    BoundaryCondition(MeshData* meshData);
    virtual ~BoundaryCondition() = default;

    virtual void fillCoeffMatrix(
        const size_t nodeID, MQBasis* RBFBasis,
        Eigen::SparseMatrix<double, Eigen::RowMajor>& spMatrix) const = 0;

    virtual void fillRhsVector(const size_t nodeID, MQBasis* RBFBasis,
                               Eigen::VectorXd& rhsVec) const = 0;

   protected:
    MeshData* meshData_;
};

#endif
