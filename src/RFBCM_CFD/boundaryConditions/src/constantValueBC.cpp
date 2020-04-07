#include "constantValueBC.hpp"

ConstantValueBC::ConstantValueBC(const double constValue)
    : constValue_(constValue){};

boundaryConditionType ConstantValueBC::type()
{
    return boundaryConditionType::constantValue;
};

// constantValueBC::addNode(const int index)
// {
//     nodeIndices_.push_back(index);
// }