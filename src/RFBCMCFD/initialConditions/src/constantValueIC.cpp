#include "constantValueIC.hpp"

ConstantValueIC::ConstantValueIC(const double constValue)
    : InitialCondition(), constValue_(constValue){};

void ConstantValueIC::fillVector(const size_t nodeID,
                                 Eigen::VectorXd& vec) const
{
    vec[nodeID] = constValue_;
}