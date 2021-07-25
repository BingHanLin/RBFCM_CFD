#include "ConstantVecValueIC.hpp"
#include "MeshData.hpp"

ConstantVecValueIC::ConstantVecValueIC(const std::array<double, 3> constValue,
                                       MeshData* meshData)
    : InitialCondition(), meshData_(meshData), constValue_(constValue){};

void ConstantVecValueIC::fillVector(const size_t nodeID,
                                    Eigen::VectorXd& vec) const
{
    const double numOfNodes = meshData_->numOfNodes();

    if (vec.size() >= numOfNodes)
    {
        vec(nodeID) = constValue_[0];
    }

    if (vec.size() >= numOfNodes * 2)
    {
        vec(numOfNodes + nodeID) = constValue_[1];
    }

    if (vec.size() >= numOfNodes * 3)
    {
        vec(numOfNodes * 2 + nodeID) = constValue_[2];
    }
}