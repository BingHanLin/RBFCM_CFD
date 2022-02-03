#ifndef CONSTANTVECVALUEBC_HPP
#define CONSTANTVECVALUEBC_HPP
#include "boundaryCondition.hpp"

class MeshData;
class MQBasis;
class ConstantVecValueBC : public BoundaryCondition
{
   public:
    ConstantVecValueBC(const std::array<double, 3> constValue,
                       MeshData* meshData);

    void fillCoeffMatrix(
        const size_t nodeID, MQBasis* RBFBasis,
        Eigen::SparseMatrix<double, Eigen::RowMajor>& spMatrix) const override;

    void fillRhsVector(const size_t nodeID, MQBasis* RBFBasis,
                       Eigen::VectorXd& rhsVec) const override;

   private:
    const std::array<double, 3> constValue_;
};

#endif