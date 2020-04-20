#include "constantValueBC.hpp"

ConstantValueBC::ConstantValueBC(const double constValue)
    : BoundaryCondition(), constValue_(constValue){};

void ConstantValueBC::fillCoeffMatrix(
    const int nodeID, const std::vector<int>& nodesCloudID,
    const std::vector<vec3d<double>>& nodesCloud,
    std::shared_ptr<MQBasis> RBFBasis, Eigen::SparseMatrix<double>& spMatrix)
{
    Eigen::VectorXd localVector =
        RBFBasis->collectOnNodes(nodesCloud, rbfOperatorType::Laplace);

    for (int i = 0; i < nodesCloudID.size(); i++)
    {
        spMatrix.insert(nodeID, nodesCloudID[i]) = localVector[i];
    }
}

boundaryConditionType ConstantValueBC::type()
{
    return boundaryConditionType::constantValue;
};

// void ConstantValueBC::setRHSValue(double &oneRHS){

// };

// constantValueBC::addNode(const int index)
// {
//     nodeIndices_.push_back(index);
// }