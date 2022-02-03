#ifndef CONSTANTVECVALUEIC_HPP
#define CONSTANTVECVALUEIC_HPP

#include "initialCondition.hpp"

class MeshData;
class ConstantVecValueIC : public InitialCondition
{
   public:
    ConstantVecValueIC(const std::array<double, 3> constValue,
                       MeshData* meshData);

    void fillVector(const size_t nodeID, Eigen::VectorXd& vec) const override;

   private:
    MeshData* meshData_;
    const std::array<double, 3> constValue_;
};

#endif
