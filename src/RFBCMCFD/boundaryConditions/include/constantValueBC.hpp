#ifndef CONSTANTVALUEBC_HPP
#define CONSTANTVALUEBC_HPP
#include "boundaryCondition.hpp"

class MeshData;
class MQBasis;
class ConstantValueBC : public BoundaryCondition
{
   public:
    ConstantValueBC(const double constValue, MeshData* meshData);

    void fillCoeffMatrix(
        const size_t nodeID, MQBasis* RBFBasis,
        Eigen::SparseMatrix<double, Eigen::RowMajor>& spMatrix) const override;

    void fillRhsVector(const size_t nodeID, MQBasis* RBFBasis,
                       Eigen::VectorXd& rhsVec) const override;

   private:
    double constValue_;
};

#endif