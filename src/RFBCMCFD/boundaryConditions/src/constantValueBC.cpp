#include "constantValueBC.hpp"
#include <iostream>

ConstantValueBC::ConstantValueBC(const double constValue, MeshData* meshData)
    : BoundaryCondition(meshData), constValue_(constValue){};

void ConstantValueBC::fillCoeffMatrix(
    const size_t nodeID, std::shared_ptr<MQBasis> RBFBasis,
    Eigen::SparseMatrix<double>& spMatrix) const
{
    spMatrix.insert(nodeID, nodeID) = 1.0;
}

void ConstantValueBC::fillRhsVector(const size_t nodeID,
                                    std::shared_ptr<MQBasis> RBFBasis,
                                    Eigen::VectorXd& rhsVec) const
{
    rhsVec(nodeID) = constValue_;
}