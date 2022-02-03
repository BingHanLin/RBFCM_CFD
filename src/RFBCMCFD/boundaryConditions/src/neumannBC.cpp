#include "NeumannBC.hpp"
#include "MQBasis.hpp"
#include "MeshData.hpp"

#include <iostream>

NeumannBC::NeumannBC(const double rhsValue, MeshData* meshData)
    : BoundaryCondition(meshData), rhsValue_(rhsValue_){};

void NeumannBC::fillCoeffMatrix(
    const size_t nodeID, MQBasis* RBFBasis,
    Eigen::SparseMatrix<double, Eigen::RowMajor>& spMatrix) const
{
    const auto cloud = meshData_->cloudByID(nodeID);

    Eigen::VectorXd localVector =
        RBFBasis->collectOnNodes(nodeID, rbfOperatorType::NEUMANN);

    for (size_t i = 0; i < localVector.size(); i++)
    {
        spMatrix.insert(nodeID, cloud.ids_[i]) = localVector[i];
    }
}

void NeumannBC::fillRhsVector(const size_t nodeID, MQBasis* RBFBasis,
                              Eigen::VectorXd& rhsVec) const
{
    rhsVec(nodeID) = rhsValue_;
}