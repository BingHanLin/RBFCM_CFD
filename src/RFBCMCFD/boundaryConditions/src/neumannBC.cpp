#include "neumannBC.hpp"
#include "MQBasis.hpp"
#include "meshData.hpp"

#include <iostream>

neumannBC::neumannBC(const double rhsValue, MeshData* meshData)
    : BoundaryCondition(meshData), rhsValue_(rhsValue_){};

void neumannBC::fillCoeffMatrix(const size_t nodeID, MQBasis* RBFBasis,
                                Eigen::SparseMatrix<double>& spMatrix) const
{
    const auto cloud = meshData_->cloudByID(nodeID);

    Eigen::VectorXd localVector =
        RBFBasis->collectOnNodes(nodeID, rbfOperatorType::NEUMANN);

    for (size_t i = 0; i < localVector.size(); i++)
    {
        spMatrix.insert(nodeID, cloud.ids_[i]) = localVector[i];
    }
}

void neumannBC::fillRhsVector(const size_t nodeID, MQBasis* RBFBasis,
                              Eigen::VectorXd& rhsVec) const
{
    rhsVec(nodeID) = rhsValue_;
}