#include "ConstantVecValueBC.hpp"
#include "MQBasis.hpp"
#include "MeshData.hpp"

#include <iostream>

ConstantVecValueBC::ConstantVecValueBC(const std::array<double, 3> constValue,
                                       MeshData* meshData)
    : BoundaryCondition(meshData), constValue_(constValue){};

void ConstantVecValueBC::fillCoeffMatrix(
    const size_t nodeID, MQBasis* RBFBasis,
    Eigen::SparseMatrix<double, Eigen::RowMajor>& spMatrix) const
{
    spMatrix.insert(nodeID, nodeID) = 1.0;
}

void ConstantVecValueBC::fillRhsVector(const size_t nodeID, MQBasis* RBFBasis,
                                       Eigen::VectorXd& rhsVec) const
{
    const double numOfNodes = meshData_->numOfNodes();

    if (rhsVec.size() >= numOfNodes)
    {
        rhsVec(nodeID) = constValue_[0];
    }

    if (rhsVec.size() >= numOfNodes * 2)
    {
        rhsVec(numOfNodes + nodeID) = constValue_[1];
    }

    if (rhsVec.size() >= numOfNodes * 3)
    {
        rhsVec(numOfNodes * 2 + nodeID) = constValue_[2];
    }
}