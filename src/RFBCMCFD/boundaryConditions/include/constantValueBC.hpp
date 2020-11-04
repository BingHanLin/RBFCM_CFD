#ifndef CONSTANTVALUEBC_HPP
#define CONSTANTVALUEBC_HPP

#include "boundaryCondition.hpp"

class meshData;

class ConstantValueBC : public BoundaryCondition
{
   public:
    ConstantValueBC(const double constValue_);
    // addNode(const size_t index);

    // boundaryConditionType type() override;

    void fillCoeffMatrix(const size_t nodeID, const nodesCloud& cloud,
                         std::shared_ptr<MQBasis> RBFBasis,
                         Eigen::SparseMatrix<double>& spMatrix) const;

    void fillRhsVector(const size_t nodeID, const nodesCloud& cloud,
                       std::shared_ptr<MQBasis> RBFBasis,
                       Eigen::VectorXd& rhsVec) const;

    // void setRHSValue(double& oneRHS);

   private:
    // std::vector<size_t> nodeIndices_;
    double constValue_;
};

#endif