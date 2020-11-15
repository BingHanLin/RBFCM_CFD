#ifndef neumannBC_HPP
#define CONSTANTVALUEBC_HPP
#include "boundaryCondition.hpp"
#include "meshData.hpp"

class neumannBC : public BoundaryCondition
{
   public:
    neumannBC(const double rhsValue, MeshData* meshData);

    void fillCoeffMatrix(const size_t nodeID, std::shared_ptr<MQBasis> RBFBasis,
                         Eigen::SparseMatrix<double>& spMatrix) const override;

    void fillRhsVector(const size_t nodeID, std::shared_ptr<MQBasis> RBFBasis,
                       Eigen::VectorXd& rhsVec) const override;

   private:
    double rhsValue_;
};

#endif