#ifndef NEUMANNBC_HPP
#define NEUMANNBC_HPP

#include "BoundaryCondition.hpp"

class MeshData;
class MQBasis;
class NeumannBC : public BoundaryCondition
{
   public:
    NeumannBC(const double rhsValue, MeshData* meshData);

    void fillCoeffMatrix(
        const size_t nodeID, MQBasis* RBFBasis,
        Eigen::SparseMatrix<double, Eigen::RowMajor>& spMatrix) const override;

    void fillRhsVector(const size_t nodeID, MQBasis* RBFBasis,
                       Eigen::VectorXd& rhsVec) const override;

   private:
    double rhsValue_;
};

#endif