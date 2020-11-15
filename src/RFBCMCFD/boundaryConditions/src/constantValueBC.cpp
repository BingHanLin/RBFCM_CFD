#include "constantValueBC.hpp"
#include <iostream>

ConstantValueBC::ConstantValueBC(const double constValue, MeshData* meshData)
    : BoundaryCondition(meshData), constValue_(constValue){};

void ConstantValueBC::fillCoeffMatrix(
    const size_t nodeID, std::shared_ptr<MQBasis> RBFBasis,
    Eigen::SparseMatrix<double>& spMatrix) const
{
    // Eigen::VectorXd localVector =
    //     RBFBasis->collectOnNodes(nodesCloud, rbfOperatorType::LAPLACE);

    // for (size_t i = 0; i < nodesCloudID.size(); i++)
    // {
    //     spMatrix.insert(nodeID, nodesCloudID[i]) = localVector[i];
    // }

    spMatrix.insert(nodeID, nodeID) = 1.0;
}

void ConstantValueBC::fillRhsVector(const size_t nodeID,
                                    std::shared_ptr<MQBasis> RBFBasis,
                                    Eigen::VectorXd& rhsVec) const
{
    rhsVec(nodeID) = constValue_;
}